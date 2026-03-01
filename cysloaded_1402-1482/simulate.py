import os
import sys

import openmm
from openff import toolkit
from openmm import app, unit

from openff.pablo import STD_CCD_CACHE, ResidueDefinition, topology_from_pdb
from openff.toolkit import ForceField

replicate = int(sys.argv[1])
print(f"This is replicate #{replicate}.")

timestep = 2 * unit.femtosecond
runtime = 10 * unit.nanoseconds
temperature = 298.15 * unit.kelvin

serine_resdef = STD_CCD_CACHE["SER"][0]
cysPPT_resdef = ResidueDefinition.anon_from_smiles("O[P@@](=O)(OCC([C@H](C(=O)NCCC(=O)NCCNC(=O)[C@H](CS)N)O)(C)C)[O-]")

adduct_resdef = ResidueDefinition.react(
    reactants=[serine_resdef, cysPPT_resdef],
    reactant_smarts=["[CH2:1][O:2][H:3]",
                     "[H:4][O:5][P:6]"],
    product_smarts=[
        "[C:1][O:2][P:6]", # Phosphorylated serine sidechain
        "[H:3][O:5][H:4]", # Water
    ],
    product_residue_names = ["PPT"],
    product_linking_bonds = [serine_resdef.linking_bond],
)[0][0]

topology = topology_from_pdb("PCP1_cysPPT.pdb",
                            additional_definitions=[adduct_resdef])

ff = ForceField(
    "openff_no_water-3.0.0-alpha0.offxml",
    "tip3p_fb.offxml",
)

interchange = ff.create_interchange(topology)

modeller = app.Modeller(interchange.topology.to_openmm(), interchange.positions.to_openmm())
modeller.addSolvent(
    ff,
    padding=1.2 * unit.nanometers,
    model="tip3p",
    ionicStrength=0.150 * unit.molar,
    boxShape="dodecahedron",
)
modeller.addExtraParticles(ff)

print("Creating system...")
system = ff.createSystem(
    modeller.topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1 * unit.nanometer,
    constraints=app.HBonds,
    rigidWater=True,
)
assert False
with open(f"{replicate}/system.xml", "w") as f:
    f.write(openmm.XmlSerializer.serialize(system))

integrator = openmm.LangevinMiddleIntegrator(temperature, 1 / unit.picosecond, timestep)
# Use the int `replicate` as the random number seed
integrator.setRandomNumberSeed(replicate)
platform = openmm.Platform.getPlatformByName("CUDA")
print("CUDA_VISIBLE_DEVICES:", os.environ.get("CUDA_VISIBLE_DEVICES"))
properties = {"Precision": "mixed"}
simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)

simulation.context.setPositions(modeller.positions)
print("Minimizing energy...")
simulation.minimizeEnergy(tolerance=5.5 * unit.kilojoules_per_mole / unit.nanometer)
positions = simulation.context.getState(positions=True).getPositions()
app.PDBFile.writeFile(
    simulation.topology, positions, open(f"{replicate}/minimized.pdb", "w")
)
simulation.saveState(f"{replicate}/minimized_state.xml")

print("Simulating production...")
simulation.reporters.append(
    app.StateDataReporter(
        f"{replicate}/production.log",
        reportInterval=int(20 * unit.picosecond / timestep),
        step=True,
        time=True,
        totalEnergy=True,
        potentialEnergy=True,
        kineticEnergy=True,
        temperature=True,
        speed=True,
    )
)
simulation.reporters.append(
    app.StateDataReporter(
        sys.stdout,
        reportInterval=int(100 * unit.picosecond / timestep),
        step=True,
        time=True,
        totalEnergy=True,
        potentialEnergy=True,
        kineticEnergy=True,
        temperature=True,
        speed=True,
        totalSteps=int(runtime / timestep),
        progress=True,
        remainingTime=True,
    )
)

simulation.reporters.append(
    app.DCDReporter(
        f"{replicate}/production.dcd",
        reportInterval=int(20 * unit.picosecond / timestep),
        enforcePeriodicBox=True,
    )
)

simulation.reporters.append(
    app.CheckpointReporter(
        f"{replicate}/checkpoint.chk", int(10 * unit.nanoseconds / timestep)
    )
)

simulation.step(int(runtime / timestep))
