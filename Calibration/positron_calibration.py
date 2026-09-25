import random
import numpy as np
import geant4_pybind as g4

OUTPUT_FILE = "junolite_positron_calibration.npz"
ENERGIES = np.arange(1, 11)
EVENTS_PER_ENERGY = 100
PMT_QE = 0.50

events = []


# =================================================================
# Detector Construction
# =================================================================
class DetectorConstruction(g4.G4VUserDetectorConstruction):

    def Construct(self):
        nist = g4.G4NistManager.Instance()
        world_mat = nist.FindOrBuildMaterial("G4_AIR")
        scintillator_mat = nist.FindOrBuildMaterial("G4_WATER")

        # Optical photon energy grid + property tables
        photon_energy = g4.G4doubleVector([2.0 * g4.eV, 2.5 * g4.eV, 3.0 * g4.eV])
        flat = lambda value: g4.G4doubleVector([value, value, value])

        world_mpt = g4.G4MaterialPropertiesTable()
        world_mpt.AddProperty("RINDEX", photon_energy, flat(1.0))
        world_mat.SetMaterialPropertiesTable(world_mpt)

        scint_mpt = g4.G4MaterialPropertiesTable()
        scint_mpt.AddProperty("RINDEX", photon_energy, flat(1.5))
        scint_mpt.AddProperty("ABSLENGTH", photon_energy, flat(10.0 * g4.m))
        scint_mpt.AddProperty("SCINTILLATIONCOMPONENT1", photon_energy, flat(1.0))
        scint_mpt.AddConstProperty("SCINTILLATIONYIELD", 1000.0 / g4.MeV)
        scint_mpt.AddConstProperty("RESOLUTIONSCALE", 1.0)
        scint_mpt.AddConstProperty("SCINTILLATIONTIMECONSTANT1", 10.0 * g4.ns)
        scintillator_mat.SetMaterialPropertiesTable(scint_mpt)

        # World
        world_size = 2.0 * g4.m
        world_solid = g4.G4Box("World", world_size / 2, world_size / 2, world_size / 2)
        world_logical = g4.G4LogicalVolume(world_solid, world_mat, "World")
        world_physical = g4.G4PVPlacement(None, g4.G4ThreeVector(0, 0, 0), world_logical,
                                           "World", None, False, 0, True)

        # Scintillator
        radius = 50.0 * g4.cm
        half_length = 50.0 * g4.cm
        scintillator_solid = g4.G4Tubs("Scintillator", 0, radius, half_length, 0, 360 * g4.deg)
        scintillator_logical = g4.G4LogicalVolume(scintillator_solid, scintillator_mat, "Scintillator")
        g4.G4PVPlacement(None, g4.G4ThreeVector(0, 0, 0), scintillator_logical,
                          "Scintillator", world_logical, False, 0, True)

        # PMT
        pmt_radius = 10.0 * g4.cm
        pmt_half_length = 1.0 * g4.cm
        pmt_solid = g4.G4Tubs("PMT", 0, pmt_radius, pmt_half_length, 0, 360 * g4.deg)
        pmt_logical = g4.G4LogicalVolume(pmt_solid, world_mat, "PMT")
        g4.G4PVPlacement(None, g4.G4ThreeVector(0, 0, 51.0 * g4.cm), pmt_logical,
                          "PMT", world_logical, False, 0, True)

        return world_physical


# =================================================================
# Primary Generator
# =================================================================
class PrimaryGeneratorAction(g4.G4VUserPrimaryGeneratorAction):

    def __init__(self):
        super().__init__()
        self.particle_gun = g4.G4ParticleGun(1)
        self.particle_gun.SetParticleByName("e+")
        self.particle_gun.SetParticlePosition(g4.G4ThreeVector(0, 0, -4.0 * g4.cm))
        self.particle_gun.SetParticleMomentumDirection(g4.G4ThreeVector(0, 0, 1))
        self.SetEnergy(1.0)

    def SetEnergy(self, energy):
        self.energy = energy
        self.particle_gun.SetParticleEnergy(energy * g4.MeV)

    def GeneratePrimaries(self, event):
        self.particle_gun.GeneratePrimaryVertex(event)


# =================================================================
# Stepping Action
# =================================================================
class SteppingAction(g4.G4UserSteppingAction):

    def __init__(self):
        super().__init__()
        self.reset()

    def reset(self):
        self.detected_pe = 0
        self.optical_photons = 0
        self.pmt_hits = 0
        self.deposited_energy = 0.0
        self.pmt_deposited_energy = 0.0

        # Positron annihilation diagnostics
        self.gamma_count = 0
        self.annihilation_gamma_count = 0
        self.annihilation_gamma_energy = 0.0
        self.gamma_deposited_energy = 0.0

    def UserSteppingAction(self, step):
        track = step.GetTrack()
        particle_name = track.GetDefinition().GetParticleName()

        pre_volume = step.GetPreStepPoint().GetTouchable().GetVolume()
        if pre_volume is not None:
            volume_name = pre_volume.GetName()
            edep = step.GetTotalEnergyDeposit()

            if volume_name == "Scintillator":
                self.deposited_energy += edep
                if particle_name == "gamma":
                    self.gamma_deposited_energy += edep
            elif volume_name == "PMT":
                self.pmt_deposited_energy += edep

        # Gamma diagnostics: count each track once, flag annihilation gammas
        if particle_name == "gamma":
            if track.GetCurrentStepNumber() == 1:
                self.gamma_count += 1
                creator = track.GetCreatorProcess()
                if creator is not None and "annihil" in creator.GetProcessName().lower():
                    self.annihilation_gamma_count += 1
                    self.annihilation_gamma_energy += track.GetKineticEnergy() / g4.MeV
            return

        if particle_name != "opticalphoton":
            return

        if track.GetCurrentStepNumber() == 1:
            self.optical_photons += 1

        post_volume = step.GetPostStepPoint().GetTouchable().GetVolume()
        if post_volume is not None and post_volume.GetName() == "PMT":
            self.pmt_hits += 1
            if random.random() < PMT_QE:
                self.detected_pe += 1
            track.SetTrackStatus(g4.fStopAndKill)


# =================================================================
# Event Action
# =================================================================
class EventAction(g4.G4UserEventAction):

    def __init__(self, stepping_action, primary_generator):
        super().__init__()
        self.sa = stepping_action
        self.pg = primary_generator

    def BeginOfEventAction(self, event):
        self.sa.reset()

    def EndOfEventAction(self, event):
        sa = self.sa

        events.append({
            "true_energy": self.pg.energy,
            "deposited_energy": sa.deposited_energy / g4.MeV,
            "pmt_deposited_energy": sa.pmt_deposited_energy / g4.MeV,
            "n_photons": sa.optical_photons,
            "pmt_hits": sa.pmt_hits,
            "n_pe": sa.detected_pe,
            "gamma_count": sa.gamma_count,
            "annihilation_gamma_count": sa.annihilation_gamma_count,
            "annihilation_gamma_energy": sa.annihilation_gamma_energy,
            "gamma_deposited_energy": sa.gamma_deposited_energy,
        })


# =================================================================
# Run
# =================================================================
run_manager = g4.G4RunManager()
run_manager.SetUserInitialization(DetectorConstruction())

physics = g4.FTFP_BERT()
physics.RegisterPhysics(g4.G4OpticalPhysics())
run_manager.SetUserInitialization(physics)

primary_generator = PrimaryGeneratorAction()
stepping_action = SteppingAction()
event_action = EventAction(stepping_action, primary_generator)

run_manager.SetUserAction(primary_generator)
run_manager.SetUserAction(stepping_action)
run_manager.SetUserAction(event_action)

run_manager.Initialize()

ui = g4.G4UImanager.GetUIpointer()
ui.ApplyCommand("/vis/open OGL")
ui.ApplyCommand("/vis/drawVolume")
ui.ApplyCommand("/vis/scene/add/trajectories")
ui.ApplyCommand("/vis/scene/endOfEventAction accumulate")

# ---------------------------------------------------------------
# Positron calibration scan
# ---------------------------------------------------------------
for energy in ENERGIES:
    print(f"\n{'-' * 40}\nRunning positron energy: {energy} MeV\n{'-' * 40}")
    primary_generator.SetEnergy(float(energy))
    run_manager.BeamOn(EVENTS_PER_ENERGY)

# ---------------------------------------------------------------
# Save results
# ---------------------------------------------------------------
data = {key: np.array([ev[key] for ev in events]) for key in events[0]}
np.savez(OUTPUT_FILE, **data)

# ---------------------------------------------------------------
# Positron annihilation diagnostic (at 10 MeV)
# ---------------------------------------------------------------
mask_10 = data["true_energy"] == 10.0

print(f"\n{'-' * 40}\nPositron Annihilation Diagnostic\n{'-' * 40}")
print("Mean gamma count:", np.mean(data["gamma_count"][mask_10]))
print("Mean annihilation gamma count:", np.mean(data["annihilation_gamma_count"][mask_10]))
print("Mean annihilation gamma energy:", np.mean(data["annihilation_gamma_energy"][mask_10]), "MeV")
print("Mean gamma deposited energy:", np.mean(data["gamma_deposited_energy"][mask_10]), "MeV")

print(f"\n{'-' * 40}\nPositron calibration run complete\n{'-' * 40}")
print("Total events:", len(events))
print(f"Saved: {OUTPUT_FILE}")
