# Conclusion

## Purpose of the project

This project was built primarily to develop and demonstrate my practical
experience with **Geant4 and detector simulation**.

The main goal was not to reproduce the JUNO detector in detail or to obtain
precision physics results. Instead, I wanted to understand how a particle
detector simulation is constructed from the ground up and to gain experience
with the different stages involved in a Monte Carlo detector simulation.

I therefore chose a JUNO-inspired liquid-scintillator detector as the physics
framework in which to develop these skills.

The simulation was built progressively, starting from particle energy
deposition and scintillation, and then adding optical photon transport, PMT
response, neutron transport, calibration, and finally neutrino-energy
reconstruction.

---

## What this project investigates

The current implementation investigates several aspects of a simplified
liquid-scintillator detector:

### 1. Scintillation and optical photon transport

The simulation tracks the conversion of deposited particle energy into
scintillation photons and transports those photons through the detector using
Geant4 optical processes.

This allowed me to work with Geant4 optical physics, material optical
properties, photon tracking, and detector boundaries.

### 2. PMT response

Optical photons reaching the PMT surfaces are converted into a simplified
photoelectron response using a PMT quantum-efficiency model.

The resulting number of photoelectrons (NPE) is used as the primary detector
observable.

This provides the connection

```text
particle energy deposition
        ↓
scintillation photons
        ↓
optical photon transport
        ↓
PMT hits
        ↓
photoelectrons
