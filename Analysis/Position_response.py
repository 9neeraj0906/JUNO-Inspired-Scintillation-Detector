import numpy as np
import matplotlib.pyplot as plt

DATA_FILE = "junolite_ibd.npz"
N_BINS = 10

# ---------------------------------------------------------------
# Load simulation data
# ---------------------------------------------------------------
data = np.load(DATA_FILE)
npe = data["prompt_npe"].astype(float)
edep = data["deposited_energy"].astype(float)
x = data["vertex_x"].astype(float)
y = data["vertex_y"].astype(float)
z = data["vertex_z"].astype(float)

radius = np.sqrt(x**2 + y**2)
valid = edep > 0.0

pe_per_mev = np.zeros_like(npe)
pe_per_mev[valid] = npe[valid] / edep[valid]


# ---------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------
def scatter_plot(xvals, yvals, xlabel, ylabel, title, filename):
    plt.figure(figsize=(8, 6))
    plt.scatter(xvals, yvals, s=12, alpha=0.5)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.close()


def binned_stats(coord, values, mask, n_bins, label):
    """Bin `values` by `coord` and print/return mean+std per bin."""
    lo, hi = np.min(coord[mask]), np.max(coord[mask])
    edges = np.linspace(lo, hi, n_bins + 1)

    centers, means, stds, counts = [], [], [], []
    for i in range(n_bins):
        bin_mask = mask & (coord >= edges[i]) & (coord < edges[i + 1])
        binned = values[bin_mask]
        if len(binned) == 0:
            continue
        centers.append(0.5 * (edges[i] + edges[i + 1]))
        means.append(np.mean(binned))
        stds.append(np.std(binned))
        counts.append(len(binned))

    centers, means, stds, counts = map(np.array, (centers, means, stds, counts))

    print(f"\nBINNED LIGHT YIELD VS {label.upper()}")
    print("-" * 60)
    print("{:>10} {:>15} {:>15} {:>10}".format(f"{label} [cm]", "Mean PE/MeV", "Std PE/MeV", "Events"))
    for c, m, s, n in zip(centers, means, stds, counts):
        print(f"{c:10.2f} {m:15.3f} {s:15.3f} {n:10d}")

    return centers, means, stds


def binned_plot(centers, means, stds, xlabel, title, filename):
    plt.figure(figsize=(8, 6))
    plt.errorbar(centers, means, yerr=stds, fmt="o", capsize=4)
    plt.xlabel(xlabel)
    plt.ylabel("Prompt NPE / Deposited Energy [PE/MeV]")
    plt.title(title)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.close()


# ---------------------------------------------------------------
# Basic statistics
# ---------------------------------------------------------------
print("=" * 60)
print("JUNO-LITE POSITION RESPONSE ANALYSIS")
print("=" * 60)
print("Total events:", len(npe))
print("Events with Edep > 0:", np.sum(valid))

print("\nDEPOSITED ENERGY")
print("-" * 60)
print("Mean :", np.mean(edep[valid]))
print("Min  :", np.min(edep[valid]))
print("Max  :", np.max(edep[valid]))

print("\nPROMPT NPE")
print("-" * 60)
print("Mean :", np.mean(npe))
print("Std  :", np.std(npe))
print("Min  :", np.min(npe))
print("Max  :", np.max(npe))

print("\nLIGHT YIELD")
print("-" * 60)
print("Mean PE/MeV :", np.mean(pe_per_mev[valid]))
print("Std  PE/MeV :", np.std(pe_per_mev[valid]))
print("Min  PE/MeV :", np.min(pe_per_mev[valid]))
print("Max  PE/MeV :", np.max(pe_per_mev[valid]))

# ---------------------------------------------------------------
# Correlations
# ---------------------------------------------------------------
corr_z = np.corrcoef(z[valid], pe_per_mev[valid])[0, 1]
corr_r = np.corrcoef(radius[valid], pe_per_mev[valid])[0, 1]

print("\nPOSITION CORRELATION")
print("-" * 60)
print("Correlation PE/MeV vs z :", corr_z)
print("Correlation PE/MeV vs r :", corr_r)

# ---------------------------------------------------------------
# Scatter plots
# ---------------------------------------------------------------
scatter_plot(edep[valid], npe[valid], "Deposited Energy [MeV]", "Prompt NPE",
             "Prompt NPE vs Deposited Energy", "position_npe_vs_edep.png")

scatter_plot(z[valid], pe_per_mev[valid], "Vertex z [cm]", "Prompt NPE / Deposited Energy [PE/MeV]",
             "Light Yield vs z Position", "light_yield_vs_z.png")

scatter_plot(radius[valid], pe_per_mev[valid], "Vertex Radius r [cm]", "Prompt NPE / Deposited Energy [PE/MeV]",
             "Light Yield vs Radial Position", "light_yield_vs_radius.png")

scatter_plot(z[valid], npe[valid], "Vertex z [cm]", "Prompt NPE",
             "Prompt NPE vs z Position", "npe_vs_z.png")

# ---------------------------------------------------------------
# Binned light yield vs z and vs radius
# ---------------------------------------------------------------
z_centers, z_means, z_stds = binned_stats(z, pe_per_mev, valid, N_BINS, "z")
binned_plot(z_centers, z_means, z_stds, "Vertex z [cm]", "Binned Light Yield vs z",
            "binned_light_yield_vs_z.png")

r_centers, r_means, r_stds = binned_stats(radius, pe_per_mev, valid, N_BINS, "r")
binned_plot(r_centers, r_means, r_stds, "Vertex Radius r [cm]", "Binned Light Yield vs Radial Position",
            "binned_light_yield_vs_radius.png")

# ---------------------------------------------------------------
# Summary
# ---------------------------------------------------------------
print("\n" + "=" * 60)
print("POSITION RESPONSE ANALYSIS COMPLETE")
print("=" * 60)
print("\nGenerated plots:")
print("  position_npe_vs_edep.png")
print("  light_yield_vs_z.png")
print("  light_yield_vs_radius.png")
print("  npe_vs_z.png")
print("  binned_light_yield_vs_z.png")
print("  binned_light_yield_vs_radius.png")
print("\nThese plots will be used to determine whether a position-dependent")
print("light-collection correction is required before the final energy-resolution study.")
print("=" * 60)
