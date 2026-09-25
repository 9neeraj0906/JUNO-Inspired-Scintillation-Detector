"""
JUNO-Lite: Geant4-Python reactor antineutrino detector simulation.

Simulates a simplified inverse beta decay (IBD) event in a
JUNO-inspired liquid scintillator detector and tracks:

    reactor antineutrino
        -> positron + neutron
        -> scintillation photons
        -> PMT photon hits
        -> photoelectron response

The simulation records energy deposition, optical-photon transport,
PMT response, and neutron-capture diagnostics for later analysis.

This is an educational detector-simulation project and is not intended
to reproduce the full JUNO detector or its precision physics model.
"""
import math
import numpy as np
import geant4_pybind as g4

# ---------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------
N_EVENTS = 1000
SEED = 12345
OUTPUT_FILE = "junolite_ibd.npz"
ENABLE_VIS = False

SCINT_YIELD_PER_MEV = 1000.0
PMT_QE = 0.50
DELAYED_BOUNDARY_NS = 1000.0

# Physics constants [MeV]
ME = 0.511
MN = 939.565
MP = 938.272
MN_MP = MN - MP
IBD_THRESHOLD = 1.806
E_NU_MAX = 10.0

# Geometry
SCINT_RADIUS = 50.0 * g4.cm
SCINT_HALF_LENGTH = 50.0 * g4.cm
PMT_RADIUS = 10.0 * g4.cm
PMT_HALF_LENGTH = 1.0 * g4.cm
N_PMT_RING = 8  # 3 rings x 8 PMTs = 24 PMTs
PMT_RING_Z = [-25.0 * g4.cm, 0.0 * g4.cm, 25.0 * g4.cm]
PMT_RING_RADIUS = 51.0 * g4.cm  # just outside the scintillator

# Volume copy numbers
VOL_WORLD, VOL_SCINT, VOL_PMT = 0, 1, 2

# Particle IDs
PDG_OPTICALPHOTON = -22
PDG_NEUTRON = 2112

rng = np.random.default_rng(SEED)
events = []


# =================================================================
# Detector Construction
# =================================================================
class DetectorConstruction(g4.G4VUserDetectorConstruction):

    def __init__(self):
        super().__init__()
        self._keep = []  # keep Python references alive

    def Construct(self):
        nist = g4.G4NistManager.Instance()
        world_mat = nist.FindOrBuildMaterial("G4_AIR")
        scint_mat = nist.FindOrBuildMaterial("G4_WATER")

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
        scint_mpt.AddConstProperty("SCINTILLATIONYIELD", SCINT_YIELD_PER_MEV / g4.MeV)
        scint_mpt.AddConstProperty("RESOLUTIONSCALE", 1.0)
        scint_mpt.AddConstProperty("SCINTILLATIONTIMECONSTANT1", 10.0 * g4.ns)
        scint_mat.SetMaterialPropertiesTable(scint_mpt)

        # World
        world_size = 2.0 * g4.m
        world_solid = g4.G4Box("World", world_size / 2, world_size / 2, world_size / 2)
        world_lv = g4.G4LogicalVolume(world_solid, world_mat, "World")
        world_pv = g4.G4PVPlacement(None, g4.G4ThreeVector(0, 0, 0), world_lv, "World",
                                     None, False, VOL_WORLD, True)

        # Scintillator
        scint_solid = g4.G4Tubs("Scintillator", 0, SCINT_RADIUS, SCINT_HALF_LENGTH, 0, 360 * g4.deg)
        scint_lv = g4.G4LogicalVolume(scint_solid, scint_mat, "Scintillator")
        scint_pv = g4.G4PVPlacement(None, g4.G4ThreeVector(0, 0, 0), scint_lv, "Scintillator",
                                     world_lv, False, VOL_SCINT, True)

        # PMTs: 3 rings of N_PMT_RING, pointing radially inward
        pmt_solid = g4.G4Tubs("PMT", 0, PMT_RADIUS, PMT_HALF_LENGTH, 0, 360 * g4.deg)
        pmt_lv = g4.G4LogicalVolume(pmt_solid, world_mat, "PMT")
        pmt_pvs = []
        pmt_id = 0

        for z in PMT_RING_Z:
            for i in range(N_PMT_RING):
                phi = 2.0 * math.pi * i / N_PMT_RING
                x = PMT_RING_RADIUS * math.cos(phi)
                y = PMT_RING_RADIUS * math.sin(phi)
                position = g4.G4ThreeVector(x, y, z)

                # rotate PMT's local +z (its axis) to point inward at azimuth phi
                rotation = g4.G4RotationMatrix()
                rotation.rotateY(-90.0 * g4.deg)
                rotation.rotateZ(phi)

                pmt_pv = g4.G4PVPlacement(rotation, position, pmt_lv, f"PMT_{pmt_id}",
                                           world_lv, False, VOL_PMT, True)
                pmt_pvs.append(pmt_pv)
                pmt_id += 1

        print(f"Created {pmt_id} PMTs")

        self._keep = [world_solid, world_lv, world_pv, scint_solid, scint_lv, scint_pv,
                      pmt_solid, pmt_lv, *pmt_pvs]
        return world_pv


# =================================================================
# IBD Neutrino Spectrum Sampler
#
# Toy reactor spectrum: flux ~ E^2 exp(-E / 1.2)
# Approximate IBD cross section: sigma ~ Ee * pe
# =================================================================
class IBDSpectrumSampler:

    def __init__(self, emin=IBD_THRESHOLD, emax=E_NU_MAX, n_grid=4000):
        e = np.linspace(emin, emax, n_grid)
        flux = e**2 * np.exp(-e / 1.2)

        e_pos = e - MN_MP
        p_pos = np.sqrt(np.clip(e_pos**2 - ME**2, 0.0, None))
        sigma = np.where(e > IBD_THRESHOLD, e_pos * p_pos, 0.0)

        w = flux * sigma
        if w.sum() <= 0:
            raise RuntimeError("IBD energy distribution has zero total weight.")

        cdf = np.concatenate(([0.0], np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(e))))
        self.cdf = cdf / cdf[-1]
        self.e = e

    def sample(self):
        return float(np.interp(rng.random(), self.cdf, self.e))


# =================================================================
# Primary Generator
# =================================================================
class PrimaryGeneratorAction(g4.G4VUserPrimaryGeneratorAction):

    def __init__(self):
        super().__init__()
        self.particle_gun = g4.G4ParticleGun(1)
        self.sampler = IBDSpectrumSampler()
        self.neutrino_energy = 0.0
        self.positron_energy = 0.0
        self.positron_kinetic_energy = 0.0
        self.neutron_energy = 0.0
        self.vertex = g4.G4ThreeVector()

    @staticmethod
    def sample_vertex():
        # uniform distribution within the scintillator cylinder
        r = math.sqrt(rng.random()) * SCINT_RADIUS
        phi = rng.uniform(0.0, 2.0 * math.pi)
        z = rng.uniform(-SCINT_HALF_LENGTH, SCINT_HALF_LENGTH)
        return g4.G4ThreeVector(r * math.cos(phi), r * math.sin(phi), z)

    @staticmethod
    def random_direction():
        cos_theta = rng.uniform(-1.0, 1.0)
        sin_theta = math.sqrt(1.0 - cos_theta**2)
        phi = rng.uniform(0.0, 2.0 * math.pi)
        return g4.G4ThreeVector(sin_theta * math.cos(phi), sin_theta * math.sin(phi), cos_theta)

    def _fire(self, name, kinetic_energy):
        gun = self.particle_gun
        gun.SetParticleByName(name)
        gun.SetParticleEnergy(kinetic_energy * g4.MeV)
        gun.SetParticlePosition(self.vertex)
        gun.SetParticleMomentumDirection(self.random_direction())
        return gun

    def GeneratePrimaries(self, event):
        self.neutrino_energy = self.sampler.sample()
        self.vertex = self.sample_vertex()

        # Positron: E_e+ (total) = E_nu - (Mn - Mp)
        self.positron_energy = self.neutrino_energy - MN_MP
        self.positron_kinetic_energy = self.positron_energy - ME
        if self.positron_kinetic_energy < 0:
            raise RuntimeError("Generated positron kinetic energy is negative.")
        self._fire("e+", self.positron_kinetic_energy).GeneratePrimaryVertex(event)

        # Neutron: toy kinetic energy, 1-50 keV
        self.neutron_energy = rng.uniform(0.001, 0.050)
        self._fire("neutron", self.neutron_energy).GeneratePrimaryVertex(event)


# =================================================================
# Stepping Action
# =================================================================
class SteppingAction(g4.G4UserSteppingAction):

    def __init__(self):
        super().__init__()
        self.reset()

    def reset(self):
        self.deposited_energy = 0.0
        self.pmt_deposited_energy = 0.0
        self.optical_photons = 0
        self.prompt_photons = 0
        self.delayed_photons = 0
        self.pmt_hits = 0
        self.prompt_hits = 0
        self.delayed_hits = 0
        self.optical_to_pmt = 0
        self.optical_to_world = 0
        self.optical_other_boundary = 0

        self.neutron_capture = False
        self.neutron_escaped = False
        self.neutron_capture_time = -1.0
        self.neutron_capture_x = 0.0
        self.neutron_capture_y = 0.0
        self.neutron_capture_z = 0.0
        self.neutron_capture_energy = -1.0

    @staticmethod
    def _copy_number(step_point):
        touchable = step_point.GetTouchable()
        return touchable.GetCopyNumber(0) if touchable is not None else -1

    def UserSteppingAction(self, step):
        track = step.GetTrack()
        pdg = track.GetDefinition().GetPDGEncoding()

        if pdg == PDG_OPTICALPHOTON:
            self._handle_optical_photon(step, track)
            return

        pre = step.GetPreStepPoint()
        pre_copy = self._copy_number(pre)
        edep = step.GetTotalEnergyDeposit()

        if edep > 0.0:
            if pre_copy == VOL_SCINT:
                self.deposited_energy += edep
            elif pre_copy == VOL_PMT:
                self.pmt_deposited_energy += edep

        if pdg == PDG_NEUTRON:
            self._handle_neutron(step, track, pre_copy)

    def _handle_optical_photon(self, step, track):
        pre = step.GetPreStepPoint()
        post = step.GetPostStepPoint()

        # count each optical photon once, at creation
        if track.GetCurrentStepNumber() == 1:
            self.optical_photons += 1
            if pre.GetGlobalTime() < DELAYED_BOUNDARY_NS * g4.ns:
                self.prompt_photons += 1
            else:
                self.delayed_photons += 1

        post_copy = self._copy_number(post)

        if post_copy == VOL_PMT:
            self.pmt_hits += 1
            if post.GetGlobalTime() < DELAYED_BOUNDARY_NS * g4.ns:
                self.prompt_hits += 1
            else:
                self.delayed_hits += 1
            self.optical_to_pmt += 1
            track.SetTrackStatus(g4.fStopAndKill)
        elif post_copy == VOL_WORLD:
            self.optical_to_world += 1
            track.SetTrackStatus(g4.fStopAndKill)
        elif post_copy != VOL_SCINT:
            self.optical_other_boundary += 1

    def _handle_neutron(self, step, track, pre_copy):
        post = step.GetPostStepPoint()

        if track.GetTrackStatus() == g4.fStopAndKill:
            proc = post.GetProcessDefinedStep()
            if proc is not None and proc.GetProcessName() == "nCapture":
                pos = post.GetPosition()
                self.neutron_capture = True
                self.neutron_capture_time = post.GetGlobalTime() / g4.ns
                self.neutron_capture_x = pos.x / g4.cm
                self.neutron_capture_y = pos.y / g4.cm
                self.neutron_capture_z = pos.z / g4.cm
                self.neutron_capture_energy = post.GetKineticEnergy() / g4.MeV
                return

        if pre_copy == VOL_SCINT and not self.neutron_capture:
            if self._copy_number(post) != VOL_SCINT:
                self.neutron_escaped = True
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
        sa, pg = self.sa, self.pg

        prompt_npe = int(rng.binomial(sa.prompt_hits, PMT_QE))
        delayed_npe = int(rng.binomial(sa.delayed_hits, PMT_QE))

        events.append({
            "neutrino_energy": pg.neutrino_energy,
            "positron_total_energy": pg.positron_energy,
            "positron_kinetic_energy": pg.positron_kinetic_energy,
            "neutron_kinetic_energy": pg.neutron_energy,
            "vertex_x": pg.vertex.x / g4.cm,
            "vertex_y": pg.vertex.y / g4.cm,
            "vertex_z": pg.vertex.z / g4.cm,

            "deposited_energy": sa.deposited_energy / g4.MeV,
            "pmt_deposited_energy": sa.pmt_deposited_energy / g4.MeV,
            "n_photons": sa.optical_photons,
            "pmt_hits": sa.pmt_hits,
            "n_pe": prompt_npe + delayed_npe,
            "prompt_npe": prompt_npe,
            "delayed_npe": delayed_npe,
            "prompt_photons": sa.prompt_photons,
            "delayed_photons": sa.delayed_photons,
            "optical_to_pmt": sa.optical_to_pmt,
            "optical_to_world": sa.optical_to_world,
            "optical_other_boundary": sa.optical_other_boundary,

            "neutron_capture": int(sa.neutron_capture),
            "neutron_escaped": int(sa.neutron_escaped),
            "neutron_capture_time": sa.neutron_capture_time,
            "neutron_capture_x": sa.neutron_capture_x,
            "neutron_capture_y": sa.neutron_capture_y,
            "neutron_capture_z": sa.neutron_capture_z,
            "neutron_capture_energy": sa.neutron_capture_energy,
        })


# =================================================================
# Run
# =================================================================
try:
    g4.G4Random.setTheSeed(SEED)
except AttributeError:
    print("Warning: could not seed the Geant4 RNG; runs may not be fully reproducible.")

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

print("\nInitializing Geant4...")
run_manager.Initialize()

ui = g4.G4UImanager.GetUIpointer()
for cmd in ("/run/verbose 0", "/event/verbose 0", "/tracking/verbose 0"):
    ui.ApplyCommand(cmd)

if ENABLE_VIS:
    for cmd in ("/vis/open OGL", "/vis/drawVolume", "/vis/scene/add/trajectories",
                "/vis/scene/endOfEventAction accumulate"):
        ui.ApplyCommand(cmd)

print("\n" + "=" * 60)
print("STARTING IBD SIMULATION")
print("=" * 60 + "\n")

run_manager.BeamOn(N_EVENTS)

# ---------------------------------------------------------------
# Save + summarize
# ---------------------------------------------------------------
data = {key: np.array([ev[key] for ev in events]) for key in events[0]} if events else {}
np.savez(OUTPUT_FILE, **data)

total = len(events)
captured = int(data["neutron_capture"].sum()) if total else 0
escaped = int(data["neutron_escaped"].sum()) if total else 0

print("\n" + "=" * 60)
print("IBD SIMULATION COMPLETE")
print("=" * 60)
print(f"Total events: {total}")

if total:
    print(f"Mean E_nu: {data['neutrino_energy'].mean():.4f} MeV")
    print(f"Mean positron kinetic energy: {data['positron_kinetic_energy'].mean():.4f} MeV")
    print(f"Mean neutron kinetic energy: {data['neutron_kinetic_energy'].mean() * 1000:.4f} keV")
    print(f"Mean deposited energy: {data['deposited_energy'].mean():.4f} MeV")
    print(f"Mean optical photons: {data['n_photons'].mean():.2f}")
    print(f"Mean photons reaching PMT: {data['optical_to_pmt'].mean():.2f}")
    print(f"Mean photons reaching world: {data['optical_to_world'].mean():.2f}")
    print(f"Mean PMT hits: {data['pmt_hits'].mean():.2f}")
    print(f"Mean prompt NPE: {data['prompt_npe'].mean():.4f}")
    print(f"Mean delayed NPE: {data['delayed_npe'].mean():.4f}")

    print("\n" + "-" * 60)
    print("NEUTRON DIAGNOSTICS")
    print("-" * 60)
    print(f"Events with neutron capture: {captured} ({100.0 * captured / total:.2f}%)")
    print(f"Events with neutron escape:  {escaped} ({100.0 * escaped / total:.2f}%)")

    t = data["neutron_capture_time"]
    mask = t >= 0
    if mask.any():
        r = np.sqrt(data["neutron_capture_x"][mask]**2
                    + data["neutron_capture_y"][mask]**2
                    + data["neutron_capture_z"][mask]**2)
        print(f"Capture time [ns]: mean {t[mask].mean():.3f}, median {np.median(t[mask]):.3f}, "
              f"min {t[mask].min():.3f}, max {t[mask].max():.3f}")
        print(f"Mean capture radius: {r.mean():.3f} cm")
    else:
        print("No neutron captures recorded.")

print(f"\nSaved: {OUTPUT_FILE}")
print("=" * 60)
