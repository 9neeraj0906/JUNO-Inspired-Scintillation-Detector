import numpy as np
import matplotlib.pyplot as plt

INPUT_FILE = "junolite_calibration.npz"

# ---------------------------------------------------------------
# Load simulation data
# ---------------------------------------------------------------
data = np.load(INPUT_FILE)
true_energy = data["true_energy"]
deposited_energy = data["deposited_energy"]
n_photons = data["n_photons"]
n_pe = data["n_pe"]

energies = np.unique(true_energy)

# ---------------------------------------------------------------
# Per-energy detector response
# ---------------------------------------------------------------
mean_edep, mean_photons, mean_npe, sigma_npe, resolution = [], [], [], [], []

print("\n" + "=" * 65)
print("JUNO-Lite Detector Calibration")
print("=" * 65)
print(f"{'E (MeV)':>8} {'Edep (MeV)':>12} {'Photons':>12} {'Mean NPE':>12} {'Resolution':>12}")
print("-" * 65)

for energy in energies:
    mask = true_energy == energy
    edep, photons, npe = deposited_energy[mask], n_photons[mask], n_pe[mask]

    mean_e, mean_p, mean_pe = np.mean(edep), np.mean(photons), np.mean(npe)
    sigma_pe = np.std(npe)
    res = sigma_pe / mean_pe

    mean_edep.append(mean_e)
    mean_photons.append(mean_p)
    mean_npe.append(mean_pe)
    sigma_npe.append(sigma_pe)
    resolution.append(res)

    print(f"{energy:8.1f} {mean_e:12.3f} {mean_p:12.1f} {mean_pe:12.2f} {res:12.4f}")

mean_edep = np.array(mean_edep)
mean_photons = np.array(mean_photons)
mean_npe = np.array(mean_npe)
sigma_npe = np.array(sigma_npe)
resolution = np.array(resolution)

# ---------------------------------------------------------------
# NPE calibration fit
# ---------------------------------------------------------------
slope, intercept = np.polyfit(mean_edep, mean_npe, 1)

print("\n" + "=" * 65)
print("Energy Calibration")
print("=" * 65)
print(f"NPE = {slope:.4f} * E + {intercept:.4f}")
print(f"Photon yield = {np.mean(mean_photons / mean_edep):.1f} photons/MeV")

# ---------------------------------------------------------------
# Plot 1: NPE calibration
# ---------------------------------------------------------------
plt.figure(figsize=(7, 5))
plt.errorbar(mean_edep, mean_npe, yerr=sigma_npe, fmt="o", capsize=4, label="Simulation")

fit_x = np.linspace(0, max(mean_edep) * 1.05, 200)
fit_y = slope * fit_x + intercept
plt.plot(fit_x, fit_y, label=f"Linear fit: NPE = {slope:.2f}E + {intercept:.2f}")

plt.xlabel("Deposited Energy (MeV)")
plt.ylabel("Mean NPE")
plt.title("JUNO-Lite Energy Calibration")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig("energy_calibration.png", dpi=300)
plt.close()

# ---------------------------------------------------------------
# Plot 2: Energy resolution
# ---------------------------------------------------------------
plt.figure(figsize=(7, 5))
plt.plot(mean_edep, resolution * 100, "o-")
plt.xlabel("Deposited Energy (MeV)")
plt.ylabel(r"Energy Resolution $\sigma_E/E$ (%)")
plt.title("JUNO-Lite Energy Resolution")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("energy_resolution.png", dpi=300)
plt.close()

# ---------------------------------------------------------------
# Plot 3: NPE distributions
# ---------------------------------------------------------------
plt.figure(figsize=(8, 5))
for energy in energies:
    mask = true_energy == energy
    plt.hist(n_pe[mask], bins=20, histtype="step", linewidth=1.5, label=f"{energy:.0f} MeV")

plt.xlabel("Number of Photoelectrons (NPE)")
plt.ylabel("Events")
plt.title("NPE Response at Different Energies")
plt.grid(alpha=0.3)
plt.legend(ncol=2)
plt.tight_layout()
plt.savefig("npe_distributions.png", dpi=300)
plt.close()

# ---------------------------------------------------------------
# 10 MeV diagnostic summary
# ---------------------------------------------------------------
mask_10 = true_energy == 10.0
npe_10 = n_pe[mask_10]

print("\n" + "=" * 65)
print("10 MeV Diagnostic")
print("=" * 65)
print(f"Mean deposited energy : {np.mean(deposited_energy[mask_10]):.3f} MeV")
print(f"Mean photons          : {np.mean(n_photons[mask_10]):.1f}")
print(f"Mean NPE              : {np.mean(npe_10):.2f}")
print(f"Sigma NPE             : {np.std(npe_10):.2f}")
print(f"Resolution            : {np.std(npe_10) / np.mean(npe_10) * 100:.2f}%")

print("\nPlots saved:")
print("  energy_calibration.png")
print("  energy_resolution.png")
print("  npe_distributions.png")
print("=" * 65)
