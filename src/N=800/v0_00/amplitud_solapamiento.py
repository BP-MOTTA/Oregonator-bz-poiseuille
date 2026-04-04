#!/usr/bin/env python3
import argparse
import os
import numpy as np
import h5py
import matplotlib.pyplot as plt


def robust_percentile_limits(x, pmin=1.0, pmax=99.0):
    lo = np.percentile(x, pmin)
    hi = np.percentile(x, pmax)
    if hi <= lo:
        hi = lo + 1e-12
    return lo, hi


def estimate_baseline(field, method="min", far_exclude_frac=0.20):
    """
    field puede ser 1D (Nx,) o 2D (M,Nx)
    """
    if method == "min":
        return float(np.min(field))

    if method == "far_field_mean":
        # idea: excluir una región alrededor del máximo y promediar el resto
        # 1D: excluye ventana central alrededor de argmax
        if field.ndim == 1:
            nx = field.size
            i0 = int(np.argmax(field))
            w = max(1, int(far_exclude_frac * nx))
            mask = np.ones(nx, dtype=bool)
            mask[max(0, i0-w):min(nx, i0+w+1)] = False
            if np.any(mask):
                return float(np.mean(field[mask]))
            return float(np.min(field))
        # 2D: usamos perfil y-promediado para detectar x_peak, luego excluimos columnas
        elif field.ndim == 2:
            prof = field.mean(axis=0)
            nx = prof.size
            i0 = int(np.argmax(prof))
            w = max(1, int(far_exclude_frac * nx))
            mask = np.ones(nx, dtype=bool)
            mask[max(0, i0-w):min(nx, i0+w+1)] = False
            if np.any(mask):
                return float(np.mean(field[:, mask]))
            return float(np.min(field))
        else:
            return float(np.min(field))

    raise ValueError(f"baseline method desconocido: {method}")


def moving_average(x, w=7):
    if w is None or w <= 1:
        return np.asarray(x)
    x = np.asarray(x, dtype=float)
    if x.size < w:
        return x
    kernel = np.ones(w) / w
    return np.convolve(x, kernel, mode="same")


def main():
    ap = argparse.ArgumentParser(description="Calcula A(t) y I_uv(t) desde run.h5.")
    ap.add_argument("h5file", help="Ruta a run.h5 (ej: v0_00/run.h5)")
    ap.add_argument("--mode", choices=["2D", "1D"], default="2D",
                    help="2D usa u[k,:,:]. 1D usa perfiles promediados en y.")
    ap.add_argument("--average_y", action="store_true",
                    help="Si mode=2D, promedia en y para trabajar 1D (U(x),V(x)).")
    ap.add_argument("--use_centering", action="store_true",
                    help="Resta baseline u0,v0 y recorta a >=0 antes del solapamiento.")
    ap.add_argument("--u0_method", choices=["min", "far_field_mean"], default="min")
    ap.add_argument("--v0_method", choices=["min", "far_field_mean"], default="min")
    ap.add_argument("--eps", type=float, default=1e-14, help="Epsilon para denom.")
    ap.add_argument("--smooth", type=int, default=0, help="Ventana moving average (0=no).")
    ap.add_argument("--outdir", default=None, help="Carpeta de salida. Default: junto al h5.")
    ap.add_argument("--prefix", default=None, help="Prefijo para archivos de salida.")
    args = ap.parse_args()

    h5path = args.h5file
    if not os.path.exists(h5path):
        raise FileNotFoundError(h5path)

    base_dir = os.path.dirname(os.path.abspath(h5path))
    outdir = args.outdir or os.path.join(base_dir, "analisis_amplitud_solapamiento")
    os.makedirs(outdir, exist_ok=True)

    prefix = args.prefix
    if prefix is None:
        prefix = os.path.basename(base_dir)  # ej: v1_08

    with h5py.File(h5path, "r") as f:
        u = f["u"]          # (nsnap, M, Nx)
        v = f["v"]
        t = f["t"][:]       # (nsnap,)

        nsnap, M, Nx = u.shape

        A_list = np.empty(nsnap, dtype=float)
        Iuv_list = np.empty(nsnap, dtype=float)

        for k in range(nsnap):
            uk = u[k, :, :]     # (M, Nx)
            vk = v[k, :, :]

            if args.mode == "1D" or args.average_y:
                work_u = uk.mean(axis=0)  # (Nx,)
                work_v = vk.mean(axis=0)
            else:
                work_u = uk               # (M,Nx)
                work_v = vk

            # A(t) = max(u)
            A_list[k] = float(np.max(work_u))

            # centering opcional
            if args.use_centering:
                u0 = estimate_baseline(work_u, method=args.u0_method)
                v0 = estimate_baseline(work_v, method=args.v0_method)
                u_eff = np.maximum(work_u - u0, 0.0)
                v_eff = np.maximum(work_v - v0, 0.0)
            else:
                u_eff = work_u
                v_eff = work_v

            num = float(np.sum(u_eff * v_eff))
            den = float(np.sum(u_eff)) + args.eps
            Iuv_list[k] = num / den

    # smoothing opcional
    A_s = moving_average(A_list, args.smooth) if args.smooth else A_list
    I_s = moving_average(Iuv_list, args.smooth) if args.smooth else Iuv_list

    # CSV
    csv_path = os.path.join(outdir, f"{prefix}_A_Iuv.csv")
    np.savetxt(
        csv_path,
        np.column_stack([t, A_list, Iuv_list]),
        delimiter=",",
        header="t,A,Iuv",
        comments=""
    )

    # Plot A(t)
    plt.figure(figsize=(9, 4))
    plt.plot(t, A_list, lw=1.5, label="A(t)=max(u)")
    if args.smooth:
        plt.plot(t, A_s, lw=2.0, label=f"MA(w={args.smooth})")
    plt.xlabel("t")
    plt.ylabel("A(t)")
    plt.title(f"{prefix} | Amplitud del activador")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.savefig(os.path.join(outdir, f"{prefix}_A.png"), dpi=150)
    plt.close()

    # Plot Iuv(t)
    plt.figure(figsize=(9, 4))
    plt.plot(t, Iuv_list, lw=1.5, label="I_uv(t)")
    if args.smooth:
        plt.plot(t, I_s, lw=2.0, label=f"MA(w={args.smooth})")
    plt.xlabel("t")
    plt.ylabel("I_uv(t)")
    plt.title(f"{prefix} | Solapamiento activador–inhibidor")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.savefig(os.path.join(outdir, f"{prefix}_Iuv.png"), dpi=150)
    plt.close()

    print(f"[OK] CSV: {csv_path}")
    print(f"[OK] Plots en: {outdir}")


if __name__ == "__main__":
    main()