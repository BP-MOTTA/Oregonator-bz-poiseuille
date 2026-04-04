#!/usr/bin/env python3
import argparse
import os
import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter


def robust_minmax(arr, pmin=1.0, pmax=99.0):
    a = np.asarray(arr)
    lo = np.percentile(a, pmin)
    hi = np.percentile(a, pmax)
    if hi <= lo:
        hi = lo + 1e-12
    return lo, hi


def attr_scalar(attrs, key, default):
    """Lee attrs HDF5 que pueden venir como escalar o array(1)."""
    if key not in attrs:
        return default
    x = np.asarray(attrs[key])
    return x.item() if x.size == 1 else float(x.flat[0])


def main():
    ap = argparse.ArgumentParser(description="Video desde run.h5 (u/v + perfil 1D en j).")
    ap.add_argument("h5file", help="Ruta a run.h5 (ej: v0_00/run.h5)")
    ap.add_argument("--field", choices=["u", "v"], default="u", help="Qué campo graficar.")
    ap.add_argument("--out", default=None, help="Salida. Default: <carpeta>/<field>.mp4")
    ap.add_argument("--fps", type=int, default=20, help="Frames por segundo.")
    ap.add_argument("--stride", type=int, default=1, help="1 de cada 'stride' snapshots.")
    ap.add_argument("--gif", action="store_true", help="Exporta GIF en vez de MP4.")
    ap.add_argument("--pmin", type=float, default=1.0, help="Percentil mínimo color.")
    ap.add_argument("--pmax", type=float, default=99.0, help="Percentil máximo color.")
    ap.add_argument("--j", type=int, default=None,
                    help="Índice j (1..M) para el perfil 1D. Default: j = M//2 (centro).")
    ap.add_argument("--figsize", nargs=2, type=float, default=[10.0, 6.0],
                    help="Tamaño figura, ej: --figsize 12 6")
    args = ap.parse_args()

    h5path = args.h5file
    if not os.path.exists(h5path):
        raise FileNotFoundError(f"No existe: {h5path}")

    out = args.out
    if out is None:
        base_dir = os.path.dirname(os.path.abspath(h5path))
        ext = "gif" if args.gif else "mp4"
        out = os.path.join(base_dir, f"{args.field}.{ext}")

    with h5py.File(h5path, "r") as f:
        if args.field not in f:
            raise KeyError(f"No existe dataset '{args.field}' en {h5path}. Keys: {list(f.keys())}")

        data = f[args.field]              # esperado: (nsnap, M, Nx)
        t = f["t"][:] if "t" in f else None

        nsnap, M, Nx = data.shape
        stride = max(1, args.stride)
        idx = np.arange(0, nsnap, stride, dtype=int)

        hx = attr_scalar(f.attrs, "hx", 1.0)
        hy = attr_scalar(f.attrs, "hy", 1.0)
        Lx = Nx * hx
        Ly = M * hy

        # j para perfil: el usuario lo da en 1..M (estilo Fortran)
        if args.j is None:
            j0 = M // 2
        else:
            if not (1 <= args.j <= M):
                raise ValueError(f"--j debe estar en [1,{M}] y diste {args.j}")
            j0 = args.j - 1  # a 0-based

        # Escala de color fija (robusta) -> evita parpadeos
        sample = data[idx, :, :]
        vmin, vmax = robust_minmax(sample, args.pmin, args.pmax)

        # Escala Y para el perfil 1D (también robusta)
        prof_sample = sample[:, j0, :]    # (nframes, Nx)
        ymin, ymax = robust_minmax(prof_sample, 1.0, 99.5)

        # ----- figura con 2 paneles -----
        fig = plt.figure(figsize=tuple(args.figsize), constrained_layout=True)
        gs = fig.add_gridspec(2, 1, height_ratios=[2.2, 1.0])

        ax0 = fig.add_subplot(gs[0, 0])
        ax1 = fig.add_subplot(gs[1, 0])

        # Heatmap (arriba)
        im = ax0.imshow(
            data[idx[0], :, :],
            origin="lower",
            aspect="auto",              # << clave para “banner”
            extent=[0.0, Lx, 0.0, Ly],
            vmin=vmin, vmax=vmax
        )
        cbar = fig.colorbar(im, ax=ax0, pad=0.01)
        cbar.set_label(args.field)

        ax0.set_ylabel("y")
        ax0.set_xlim(0.0, Lx)
        ax0.set_ylim(0.0, Ly)

        # Perfil 1D (abajo) en j0
        x = np.linspace(0.0, Lx, Nx, endpoint=False)
        line, = ax1.plot(x, data[idx[0], j0, :], lw=2)
        ax1.set_xlabel("x")
        ax1.set_ylabel(f"{args.field}(x) @ j={j0+1}")
        ax1.set_xlim(0.0, Lx)
        ax1.set_ylim(0, 1)

        supt = fig.suptitle("")

        def update(frame_i):
            k = idx[frame_i]
            frame = data[k, :, :]
            im.set_data(frame)

            prof = data[k, j0, :]
            line.set_ydata(prof)

            if t is not None:
                supt.set_text(f"{args.field} | t={t[k]:.6f} | snap {k}/{nsnap-1} | j={j0+1}")
            else:
                supt.set_text(f"{args.field} | snap {k}/{nsnap-1} | j={j0+1}")

            return (im, line, supt)

        anim = FuncAnimation(fig, update, frames=len(idx), interval=1000.0/args.fps, blit=False)

        if args.gif:
            writer = PillowWriter(fps=args.fps)
            anim.save(out, writer=writer)
        else:
            writer = FFMpegWriter(fps=args.fps, bitrate=2000)
            anim.save(out, writer=writer)

        plt.close(fig)

    print(f"OK -> video guardado en: {out}")


if __name__ == "__main__":
    main()