#!/usr/bin/env python3
import argparse
import os
import numpy as np
import h5py
import matplotlib.pyplot as plt


def attr_scalar(attrs, key, default):
    """Lee attrs HDF5 que pueden venir como escalar o array(1)."""
    if key not in attrs:
        return default
    x = np.asarray(attrs[key])
    return x.item() if x.size == 1 else float(x.flat[0])


def circular_distance_idx(i, i0, Nx):
    d = np.abs(i - i0)
    return np.minimum(d, Nx - d)


def estimate_baseline_constant(U_ref, V_ref, strategy="constant_min"):
    """
    U_ref, V_ref: 1D (Nx,) ya promediados en y.
    """
    if strategy == "constant_min":
        return float(U_ref.min()), float(V_ref.min())

    if strategy == "constant_far_field":
        Nx = U_ref.size
        i_peak = int(np.argmax(U_ref))
        dist = circular_distance_idx(np.arange(Nx), i_peak, Nx)
        far_mask = dist > (Nx // 4)
        if np.any(far_mask):
            return float(U_ref[far_mask].mean()), float(V_ref[far_mask].mean())
        return float(U_ref.min()), float(V_ref.min())

    raise ValueError(f"baseline_strategy inválida: {strategy}")


def moving_average_reflect(x, w):
    """
    Moving average con padding reflect para evitar artefactos de borde (líneas verticales).
    Devuelve array del mismo largo que x.
    """
    x = np.asarray(x, dtype=float)
    if w <= 1:
        return x
    if x.size < w:
        return x

    h = w // 2
    kernel = np.ones(w, dtype=float) / float(w)

    # reflect evita el "zero padding" implícito que genera extremos raros
    xp = np.pad(x, pad_width=(h, h), mode="reflect")
    y = np.convolve(xp, kernel, mode="valid")  # len(y) == len(x)
    return y


def unwrap_periodic_x(x, Lx):
    """Unwrap por saltos periódicos usando umbral Lx/2."""
    x = np.asarray(x, dtype=float)
    xu = x.copy()
    off = 0.0
    for k in range(1, xu.size):
        dx = xu[k] - xu[k - 1]
        if dx > 0.5 * Lx:
            off -= Lx
        elif dx < -0.5 * Lx:
            off += Lx
        xu[k] += off
    return xu


def compute_plot_slice(n, smooth_w):
    """
    Para evitar líneas verticales por bordes:
      - siempre quita primer y último punto
      - si smooth_w>1, quita además h=w//2 puntos a cada lado
    """
    if n <= 2:
        return slice(None)

    start = 1
    end = n - 1

    if smooth_w and smooth_w > 1:
        h = smooth_w // 2
        start = max(start, h)
        end = min(end, n - h)

    if end <= start:
        return slice(None)

    return slice(start, end)


def main():
    ap = argparse.ArgumentParser(
        description="Prueba A: Iuv_local(t) centrado en xcm(t) con periodicidad en x."
    )
    ap.add_argument("h5file", help="Ruta a run.h5 (ej: v1_08/run.h5)")

    ap.add_argument("--W", type=float, required=True,
                    help="Ancho físico de ventana local W_win_phys (misma unidad que hx).")
    ap.add_argument("--baseline_strategy", choices=["constant_min", "constant_far_field"],
                    default="constant_min")
    ap.add_argument("--baseline_k", type=int, default=1,
                    help="Número de snapshots iniciales para baseline constante (promedio).")

    ap.add_argument("--use_centering", action="store_true",
                    help="Usar u_eff=max(U-u0,0), v_eff=max(V-v0,0).")
    ap.add_argument("--eps", type=float, default=1e-12)
    ap.add_argument("--smooth", type=int, default=0,
                    help="Moving average window (0=no).")
    ap.add_argument("--outdir", default=None,
                    help="Salida. Default: <carpeta>/analisis_solapamiento_local")
    ap.add_argument("--prefix", default=None,
                    help="Prefijo de archivos. Default: nombre carpeta (vX_YY).")

    args = ap.parse_args()

    h5path = args.h5file
    if not os.path.exists(h5path):
        raise FileNotFoundError(h5path)

    base_dir = os.path.dirname(os.path.abspath(h5path))
    outdir = args.outdir or os.path.join(base_dir, "analisis_solapamiento_local")
    os.makedirs(outdir, exist_ok=True)
    prefix = args.prefix or os.path.basename(base_dir)

    with h5py.File(h5path, "r") as f:
        u = f["u"]          # (nsnap, M, Nx) en tu convención actual
        v = f["v"]
        t = f["t"][:]       # (nsnap,)
        nsnap, M, Nx = u.shape

        hx = attr_scalar(f.attrs, "hx", 1.0)
        Lx = Nx * hx

        # trig para CM circular
        theta = 2.0 * np.pi * np.arange(Nx) / Nx
        cosT = np.cos(theta)
        sinT = np.sin(theta)

        # baseline constante con k snapshots iniciales (promedio)
        kref = max(1, min(args.baseline_k, nsnap))
        # promedio en y -> 1D; luego promedio en kref snapshots
        U_ref = np.mean([u[k, :, :].mean(axis=0) for k in range(kref)], axis=0)
        V_ref = np.mean([v[k, :, :].mean(axis=0) for k in range(kref)], axis=0)
        u0_const, v0_const = estimate_baseline_constant(U_ref, V_ref, args.baseline_strategy)

        # ventana local
        W_idx_raw = int(np.rint(args.W / hx))
        W_idx = max(0, min(W_idx_raw, Nx // 2 - 1))
        local_idx = np.r_[0:W_idx + 1, Nx - W_idx:Nx] if W_idx > 0 else np.array([0], dtype=int)

        xcm = np.full(nsnap, np.nan, dtype=float)
        xcm_un = np.full(nsnap, np.nan, dtype=float)
        Iglob = np.full(nsnap, np.nan, dtype=float)
        Iloc  = np.full(nsnap, np.nan, dtype=float)
        A     = np.full(nsnap, np.nan, dtype=float)

        for k in range(nsnap):
            # perfiles 1D por promedio en y (Prueba A)
            U = u[k, :, :].mean(axis=0)   # (Nx,)
            V = v[k, :, :].mean(axis=0)

            A[k] = float(U.max())

            if args.use_centering:
                u_eff = np.maximum(U - u0_const, 0.0)
                v_eff = np.maximum(V - v0_const, 0.0)
            else:
                u_eff = U
                v_eff = V

            denom = float(u_eff.sum())
            if denom <= args.eps:
                # pulso inexistente / no interpretable
                continue
            denom += args.eps

            # CM circular
            C = float(np.sum(u_eff * cosT)) / denom
            S = float(np.sum(u_eff * sinT)) / denom
            th = np.arctan2(S, C)
            if th < 0:
                th += 2.0 * np.pi
            x_cm = (Lx / (2.0 * np.pi)) * th
            xcm[k] = x_cm

            # solapamiento global
            Iglob[k] = float(np.sum(u_eff * v_eff)) / denom

            # shift circular para centrar pulso en i=0
            i_shift = int(np.rint(x_cm / hx)) % Nx
            Uc = np.roll(u_eff, -i_shift)
            Vc = np.roll(v_eff, -i_shift)

            num_loc = float(np.sum(Uc[local_idx] * Vc[local_idx]))
            den_loc = float(np.sum(Uc[local_idx])) + args.eps
            Iloc[k] = num_loc / den_loc

        # unwrap xcm (en [0,Lx)) -> continuo
        # si hay NaNs por denom pequeño, limpiamos antes
        valid = np.isfinite(xcm)
        if np.any(valid):
            xcm_un_valid = unwrap_periodic_x(xcm[valid], Lx)
            xcm_un[valid] = xcm_un_valid

    # smoothing opcional (solo para plot). Usamos reflect SIEMPRE que smooth>1
    if args.smooth and args.smooth > 1:
        Iglob_s = moving_average_reflect(Iglob, args.smooth)
        Iloc_s  = moving_average_reflect(Iloc,  args.smooth)
        A_s     = moving_average_reflect(A,     args.smooth)
        xcm_s   = moving_average_reflect(xcm_un, args.smooth)
    else:
        Iglob_s, Iloc_s, A_s, xcm_s = Iglob, Iloc, A, xcm_un

    # índice de plot que evita extremos raros (siempre quita 1er/último + bordes del smooth)
    idx_plot = compute_plot_slice(len(t), args.smooth)

    # máscara finita para no dibujar líneas raras por NaNs
    def plot_masked(ax, tt, yy, **kwargs):
        m = np.isfinite(tt) & np.isfinite(yy)
        ax.plot(tt[m], yy[m], **kwargs)

    # CSV
    csv_path = os.path.join(outdir, f"{prefix}_PruebaA.csv")
    np.savetxt(
        csv_path,
        np.column_stack([t, xcm, xcm_un, A, Iglob, Iloc]),
        delimiter=",",
        header="t,xcm,xcm_unwrapped,A,Iuv_global,Iuv_local",
        comments=""
    )

    # Plots
    fig, ax = plt.subplots(figsize=(10, 5))
    plot_masked(ax, t[idx_plot], Iglob_s[idx_plot], lw=1.8, label="I_uv_global")
    plot_masked(ax, t[idx_plot], Iloc_s[idx_plot],  lw=2.2, label=f"I_uv_local (W={args.W:g})")
    ax.set_xlabel("t"); ax.set_ylabel("I_uv")
    ax.set_title(f"{prefix} | Solapamiento global vs local (centrado en x_cm)")
    ax.grid(True, alpha=0.3); ax.legend()
    fig.savefig(os.path.join(outdir, f"{prefix}_Iuv_global_vs_local.png"), dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 4))
    plot_masked(ax, t[idx_plot], xcm_s[idx_plot], lw=2.0)
    ax.set_xlabel("t"); ax.set_ylabel("x_cm (unwrapped)")
    ax.set_title(f"{prefix} | x_cm(t) (media circular + unwrap)")
    ax.grid(True, alpha=0.3)
    fig.savefig(os.path.join(outdir, f"{prefix}_xcm_unwrapped.png"), dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 4))
    plot_masked(ax, t[idx_plot], A_s[idx_plot], lw=2.0)
    ax.set_xlabel("t"); ax.set_ylabel("A(t)=max(U)")
    ax.set_title(f"{prefix} | Amplitud A(t)")
    ax.grid(True, alpha=0.3)
    fig.savefig(os.path.join(outdir, f"{prefix}_A.png"), dpi=150)
    plt.close(fig)

    print(f"[OK] CSV  -> {csv_path}")
    print(f"[OK] Plots-> {outdir}")
    print(f"[INFO] Nx={Nx}, hx={hx}, Lx={Lx}")
    print(f"[INFO] W_idx_raw={W_idx_raw}, W_idx_used={W_idx}, smooth={args.smooth}, idx_plot={idx_plot}")


if __name__ == "__main__":
    main()