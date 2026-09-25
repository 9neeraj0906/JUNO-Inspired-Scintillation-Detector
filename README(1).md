# JUNO-Inspired Scintillation Detector

**Repository:** [github.com/9neeraj0906/JUNO-Inspired-Scintillation-Detector](https://github.com/9neeraj0906/JUNO-Inspired-Scintillation-Detector)

## Geant4-Python Monte Carlo Simulation of Reactor Antineutrino Detection

This project presents a simplified Monte Carlo simulation of a liquid-scintillator
reactor antineutrino detector developed using **Geant4** and **Python**.

The objective is to study the detector-response chain for inverse beta decay (IBD),

$$
\bar{\nu}_e + p \rightarrow e^+ + n,
$$

with particular emphasis on:

- scintillation photon production,
- optical photon transport,
- photomultiplier-tube (PMT) response,
- positron energy deposition,
- neutron transport and capture,
- delayed signals,
- detector position dependence, and
- neutrino-energy reconstruction.

The project is intended as a **Monte Carlo detector-simulation study** rather than
a precision reproduction of the JUNO experiment. The detector geometry, optical
properties, event generation, and detector response are deliberately simplified
so that the individual components of the simulation can be studied explicitly.

---

## 1. Motivation

Liquid-scintillator detectors provide a powerful method for detecting reactor
electron antineutrinos through inverse beta decay.

The IBD interaction produces a positron and a neutron. The positron deposits its
energy almost immediately and produces a prompt scintillation signal, while the
neutron undergoes thermalization followed by capture, producing a delayed signal.

The resulting prompt-delayed coincidence provides a characteristic signature
of reactor antineutrino interactions.

The detector response can therefore be represented schematically as

$$
\bar{\nu}_e
\rightarrow
e^+ + n
\rightarrow
\text{energy deposition}
\rightarrow
\text{scintillation photons}
\rightarrow
\text{optical transport}
\rightarrow
\text{PMT hits}
\rightarrow
\text{photoelectrons}.
$$

This project implements this chain in a simplified Geant4 environment and uses
the resulting detector observables to investigate how the deposited energy is
translated into a measurable photoelectron signal.

---

## 2. Detector Model

The simulated detector consists of a cylindrical liquid-scintillator volume
surrounded by an array of simplified PMTs.

The present geometry contains:

- cylindrical scintillator volume,
- 24 PMTs arranged around the detector,
- optical-photon transport,
- scintillation photon production,
- simplified optical properties,
- neutron transport,
- neutron capture,
- event-level detector-response recording.

The detector is intentionally much smaller and simpler than JUNO.

It should therefore be interpreted as a **JUNO-inspired detector model**, rather
than a geometrical or optical model of the actual JUNO detector.

---

## 3. Simulation Framework

The simulation is implemented using:

- **Geant4 11.4.2**
- **geant4_pybind 0.1.3**
- **Python 3.12**
- NumPy
- SciPy
- Matplotlib

Geant4 is responsible for particle transport and detector interactions,
including electromagnetic processes, neutron transport, and optical-photon
propagation.

Python is used to construct the simulation, generate events, collect detector
observables, and perform the subsequent analysis.

---

## 4. Event Generation

Reactor antineutrino events are generated using a simplified neutrino-energy
spectrum.

For each event, a neutrino energy is sampled and an IBD final state is
constructed consisting of:

$$
e^+ + n.
$$

The positron and neutron are generated from a common interaction vertex inside
the cylindrical detector.

The positron provides the primary prompt energy deposition, while the neutron
is transported independently through the detector until capture or escape.

The event generator is intentionally simplified and does not attempt to model
the complete reactor neutrino production spectrum or the full IBD differential
cross section.

---

## 5. Optical Photon Simulation

Energy deposited in the scintillator produces scintillation photons.

These photons are propagated through the detector using Geant4's optical
photon processes.

The simulation records:

- number of scintillation photons produced,
- photons reaching the PMT surfaces,
- photons reaching the world boundary,
- PMT photon hits,
- photoelectron response.

A simplified PMT quantum efficiency is applied to photons reaching the PMT
surfaces to obtain the number of detected photoelectrons (NPE).

Thus the simulation connects deposited particle energy to an experimentally
measurable detector observable:

$$
E_{\mathrm{dep}}
\rightarrow
N_{\gamma}
\rightarrow
N_{\mathrm{PMT}}
\rightarrow
N_{\mathrm{PE}}.
$$

---

## 6. Energy Calibration

Before reconstructing neutrino energies, the detector response is studied
using monoenergetic charged particles.

Electron and positron calibration runs are used to determine the relationship
between deposited energy and detected photoelectrons.

The detector response is approximated by

$$
N_{\mathrm{PE}} = aE_{\mathrm{dep}} + b.
$$

The calibration study provides:

- mean deposited energy,
- mean scintillation-photon yield,
- mean NPE,
- NPE fluctuations,
- energy resolution.

The resulting response is then used as the basis for neutrino-energy
reconstruction.

---

## 7. Energy Resolution

The statistical spread of the NPE response provides a measure of detector
energy resolution:

$$
\frac{\sigma_E}{E}
\approx
\frac{\sigma_{\mathrm{NPE}}}{\mathrm{NPE}}.
$$

The energy dependence of this quantity is studied using the calibration sample.

This provides a direct demonstration of how fluctuations in the detected
photoelectron signal propagate into the reconstructed energy.

---

## 8. Neutrino Energy Reconstruction

The calibrated detector response is used to reconstruct the visible energy
associated with each simulated IBD event.

The detector response is approximated by

```math
N_{\mathrm{PE}} = aE_{\mathrm{vis}} + b
```

where \(N_{\mathrm{PE}}\) is the detected number of photoelectrons and
\(E_{\mathrm{vis}}\) is the visible energy deposited in the detector.

The visible energy is reconstructed as

```math
E_{\mathrm{vis}}^{\mathrm{rec}}
=
\frac{N_{\mathrm{PE}}-b}{a}.
```
A simplified relation between visible energy and neutrino energy is then used:

$$
E_{\nu}^{\mathrm{rec}}
\approx
E_{\mathrm{vis}}^{\mathrm{rec}} + 0.78~\mathrm{MeV}.
$$

The reconstructed neutrino energy is compared directly with the generator-level
neutrino energy.

The analysis investigates:

- true versus reconstructed energy,
- reconstruction bias,
- residual distributions,
- energy resolution,
- reconstructed neutrino spectrum.

---

## 9. Position Dependence

A realistic scintillator detector does not necessarily have an identical
response at every position.

The simulation therefore studies the light yield as a function of the
interaction vertex.

The quantity

$$
\frac{N_{\mathrm{PE}}}{E_{\mathrm{dep}}}
$$

is evaluated as a function of:

- cylindrical radius $r$,
- detector coordinate $z$.

This provides a simple study of detector-response uniformity and demonstrates
how optical geometry and PMT coverage influence the measured signal.

The current detector contains PMTs primarily around the cylindrical surface,
so position-dependent variations are expected, particularly along the detector
axis.

---

## 10. Neutron Detection

The neutron produced in the IBD interaction is transported through the
detector.

The simulation records whether the neutron:

- is captured inside the detector, or
- escapes the detector volume.

For captured neutrons, the simulation records quantities such as:

- capture time,
- capture position,
- capture coordinates,
- delayed detector response.

The capture-time distribution provides a simple demonstration of the delayed
component of the IBD signature.

The delayed NPE distribution is also studied to connect neutron transport and
capture to the detector's delayed signal.

---

## 11. Results

The current simulation produces several detector-level observables.

### Energy calibration

The calibration demonstrates an approximately linear relationship between
deposited energy and detected photoelectrons over the simulated energy range.

For the present detector configuration, the electron calibration produced a
response of approximately

$$
N_{\mathrm{PE}}
\approx
22.8E_{\mathrm{dep}} + 8.3.
$$

The corresponding photoelectron fluctuations decrease in relative size as
the deposited energy increases.

These values describe this particular simplified detector configuration and
should not be interpreted as JUNO detector performance parameters.

### Neutrino energy reconstruction

For the 1000-event IBD sample, the simulation produced a mean generated
neutrino energy of approximately

$$
\langle E_\nu\rangle \approx 6.06~\mathrm{MeV}.
$$

Using the simplified detector calibration, the mean reconstructed energy was
approximately

$$
\langle E_\nu^{\mathrm{rec}}\rangle \approx 5.72~\mathrm{MeV}.
$$

The difference illustrates the effect of detector-response fluctuations,
partial energy containment, and the simplified optical/PMT model.

The reconstructed-energy distribution also contains events with significant
deviations from the generator energy. These events are associated with the
strong position dependence of the simplified detector response, particularly
near the detector boundary.

### Neutron response

The simulation also demonstrates neutron transport and delayed capture.

In the 1000-event sample:

- 66 events contained a recorded neutron capture,
- 107 events recorded neutron escape,
- the mean capture time for captured neutrons was approximately 16.3 μs.

These quantities are properties of the present simplified geometry and physics
configuration rather than predictions for JUNO.

---

## 12. Limitations

This project deliberately makes several approximations.

### Detector geometry

The detector is a small cylindrical model with a simplified PMT arrangement.
It does not reproduce the full JUNO detector geometry, PMT population, or
optical coverage.

### Scintillator properties

The present simulation uses simplified optical properties rather than a
complete material model of JUNO's liquid scintillator.

### Reactor spectrum and IBD generation

The neutrino spectrum and IBD final-state generation are simplified. The
simulation does not implement a complete reactor flux model or a precision
IBD cross-section treatment.

### PMT response

The PMT response is represented using a simplified quantum-efficiency model.
The simulation does not reproduce the detailed response of real JUNO PMTs,
electronics, dark noise, afterpulsing, or charge response.

### Energy reconstruction

The reconstruction currently uses a first-order global calibration. It does
not perform a full position-dependent detector calibration.

### Neutron capture

The neutron analysis demonstrates transport, capture, and timing, but the
present configuration is not intended to reproduce the complete neutron
capture physics and detector response of JUNO.

These limitations are intentional for the current stage of the project. They
allow the simulation and analysis chain to remain transparent while providing
a framework that can be extended toward more realistic detector modelling.

---

## 13. Project Structure

```text
JUNO-Inspired-Scintillation-Detector/
│
├── Simulation/
│   └── ibd.py
│
├── Calibration/
│   ├── electron_calibration.py
│   └── positron_calibration.py
│
├── Analysis/
│   ├── Position_response.py
│   └── neutron_analysis.py
│
├── Requirements.txt
├── LICENSE
└── README.md
```

---

## 14. Future Development

The current simulation provides a foundation for increasing the realism of
the detector model.

Possible extensions include:

- more realistic liquid-scintillator optical properties,
- wavelength-dependent refractive index and absorption,
- realistic PMT geometry and quantum efficiency,
- improved optical coverage,
- position-dependent energy calibration,
- more realistic reactor antineutrino spectra,
- complete IBD kinematics,
- improved neutron-capture modelling,
- detector response matrices,
- systematic uncertainty studies,
- neutrino oscillation probabilities,
- comparison with published reactor-neutrino spectra.

The project is therefore structured as a progressively extensible detector
simulation rather than a fixed reproduction of an existing experiment.
