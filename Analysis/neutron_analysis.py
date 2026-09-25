import numpy as np
import matplotlib.pyplot as plt

DATA_FILE = "junolite_ibd.npz"
ANALYSIS_OUTPUT = "neutron_analysis.npz"

# ---------------------------------------------------------------
# Load data
# ---------------------------------------------------------------
data = np.load(DATA_FILE)

neutron_capture = data["neutron_capture"].astype(bool)
neutron_escaped = data["neutron_escaped"].astype(bool)
capture_time = data["neutron_capture_time"].astype(float)
capture_x = data["neutron_capture_x"].astype(float)
capture_y = data["neutron_capture_y"].astype(float)
capture_z = data["neutron_capture_z"].astype(float)
capture_energy = data["neutron_capture_energy"].astype(float)
delayed_npe = data["delayed_npe"].astype(float)

n_events = len(neutron_capture)
n_capture = int(np.sum(neutron_capture))
n_escape = int(np.sum(neutron_escaped))

capture_times = capture_time[neutron_capture]
capture_x_valid = capture_x[neutron_capture]
capture_y_valid = capture_y[neutron_capture]
capture_z_valid = capture_z[neutron_capture]
capture_energy_valid = capture_energy[neutron_capture]
capture_radius = np.sqrt(capture_x_valid**2 + capture_y_valid**2 + capture_z_valid**2)


# ---------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------
def print_stats(title, values, unit="", fmt=".3f", extra=None):
    """Print a labeled block of mean/min/max (+ optional extra stats)."""
    print(f"\n{title}")
    print("-" * 60)
    if len(values) == 0:
        print("No data available.")
        return

    suffix = f" {unit}" if unit else ""
    print(f"Mean : {np.mean(values):{fmt}}{suffix}")
    if extra == "median":
        print(f"Median : {np.median(values):{fmt}}{suffix}")
    if extra == "std":
        print(f"Std  : {np.std(values):{fmt}}{suffix}")
    print(f"Min  : {np.min(values):{fmt}}{suffix}")
    print(f"Max  : {np.max(values):{fmt}}{suffix}")


def hist_plot(values, bins, xlabel, title, filename, xscale=None):
    plt.figure(figsize=(8, 6))
    plt.hist(values, bins=bins)
    if xscale:
        plt.xscale(xscale)
    plt.xlabel(xlabel)
    plt.ylabel("Events")
    plt.title(title)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.close()


# ---------------------------------------------------------------
# Summary
# ---------------------------------------------------------------
print("=" * 60)
print("JUNO-LITE NEUTRON ANALYSIS")
print("=" * 60)

print("\nEVENT SUMMARY")
print("-" * 60)
print("Total events        :", n_events)
print("Neutron captures    :", n_capture)
print("Neutron escapes     :", n_escape)
if n_events > 0:
    print(f"Capture fraction     : {100.0 * n_capture / n_events:.2f}%")
    print(f"Escape fraction      : {100.0 * n_escape / n_events:.2f}%")

print_stats("NEUTRON CAPTURE TIME", capture_times, unit="ns", extra="median")
print_stats("CAPTURE ENERGY", capture_energy_valid, unit="MeV", fmt=".6f")
print_stats("DELAYED PHOTOELECTRONS", delayed_npe, extra="std")

print("\nNEUTRON CAPTURE POSITION")
print("-" * 60)
if n_capture > 0:
    print(f"Mean x : {np.mean(capture_x_valid):.3f} cm")
    print(f"Mean y : {np.mean(capture_y_valid):.3f} cm")
    print(f"Mean z : {np.mean(capture_z_valid):.3f} cm")
    print(f"Mean r : {np.mean(capture_radius):.3f} cm")
    print(f"Min r  : {np.min(capture_radius):.3f} cm")
    print(f"Max r  : {np.max(capture_radius):.3f} cm")
else:
    print("No capture positions available.")

# ---------------------------------------------------------------
# Plots
# ---------------------------------------------------------------
if n_capture > 0:
    positive_times = capture_times[capture_times > 0]
    if len(positive_times) > 0:
        # log bins: capture times span a wide dynamic range
        bins = np.logspace(np.log10(positive_times.min()), np.log10(positive_times.max()), 25)
        hist_plot(positive_times, bins, "Neutron Capture Time [ns]",
                  "Neutron Capture-Time Distribution", "neutron_capture_time.png", xscale="log")
    else:
        hist_plot(capture_times, 20, "Neutron Capture Time [ns]",
                  "Neutron Capture-Time Distribution", "neutron_capture_time.png")

    hist_plot(capture_radius, 20, "Capture Radius [cm]",
              "Neutron Capture Position: Radius", "neutron_capture_radius.png")

    hist_plot(capture_z_valid, 20, "Capture z Position [cm]",
              "Neutron Capture Position: z", "neutron_capture_z.png")

    plt.figure(figsize=(8, 6))
    plt.scatter(capture_x_valid, capture_y_valid, s=15, alpha=0.7)
    plt.xlabel("Capture x [cm]")
    plt.ylabel("Capture y [cm]")
    plt.title("Neutron Capture Positions")
    plt.axis("equal")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("neutron_capture_xy.png", dpi=150)
    plt.close()

hist_plot(delayed_npe, 25, "Delayed NPE", "Delayed Photoelectron Distribution",
          "delayed_npe_distribution.png")

# ---------------------------------------------------------------
# Save analysis data
# ---------------------------------------------------------------
np.savez(
    ANALYSIS_OUTPUT,
    capture_time=capture_times,
    capture_x=capture_x_valid,
    capture_y=capture_y_valid,
    capture_z=capture_z_valid,
    capture_radius=capture_radius,
    capture_energy=capture_energy_valid,
    delayed_npe=delayed_npe,
)

# ---------------------------------------------------------------
# Final message
# ---------------------------------------------------------------
print("\n" + "=" * 60)
print("NEUTRON ANALYSIS COMPLETE")
print("=" * 60)

print("\nGenerated plots:")
if n_capture > 0:
    print("  neutron_capture_time.png")
    print("  neutron_capture_radius.png")
    print("  neutron_capture_z.png")
    print("  neutron_capture_xy.png")
print("  delayed_npe_distribution.png")

print(f"\nSaved:\n  {ANALYSIS_OUTPUT}")
print("=" * 60)
