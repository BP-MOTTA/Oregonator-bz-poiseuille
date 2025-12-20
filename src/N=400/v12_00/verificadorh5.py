import numpy as np
import h5py
import matplotlib.pyplot as plt

fname = "run.h5"

with h5py.File(fname, "r") as f:
    u = f["u"]          # (nsnap, M, Nx_phys)
    v = f["v"]
    t = f["t"][:]       # (nsnap,)

    print("u shape:", u.shape)
    print("v shape:", v.shape)
    print("t shape:", t.shape)

    assert u.shape == v.shape
    assert u.shape[0] == t.shape[0], "nsnap de u/v no coincide con len(t)"

    # Chequeo dt entre snapshots (debería ser ~Nsnap*dt)
    dt_snap = np.diff(t)
    print("t[0], t[-1]:", t[0], t[-1])
    print("dt_snap min/max:", dt_snap.min(), dt_snap.max())

    # valores finitos
    assert np.isfinite(u[0]).all() and np.isfinite(u[-1]).all()

    # Plot de un snapshot
    k = 1700
    plt.figure()
    plt.title(f"u snapshot {k}, t={t[k]:.6f}")
    plt.imshow(u[k, :, :], origin="lower", aspect="auto")
    plt.colorbar()
    plt.xlabel("x index (Nx_phys)")
    plt.ylabel("y index (M)")
    plt.show()