import sys
import os

import openmm
from openmm import app, unit

replicate = int(sys.argv[1])
print(f"This is replicate #{replicate}.")
#os.makedirs(str(replicate), exist_ok=True)

timestep = 2 * unit.femtosecond
runtime = 10 * unit.nanoseconds
temperature = 298.15 * unit.kelvin

pdb = app.PDBFile("5u3h_prepped.pdb")
#ff = app.ForceField("amber19/protein.ff19SB.xml", "amber19/tip3pfb.xml")
ff = app.ForceField("protein.ff19SB.xml", "tip3pfb.xml")
modeller = app.Modeller(pdb.topology, pdb.positions)
modeller.addSolvent(
    ff, padding=1.2 * unit.nanometers, model="tip3p", ionicStrength=0.150 * unit.molar, boxShape="dodecahedron"
)
modeller.addExtraParticles(ff)
# print(modelddler.topology.getPeriodicBoxVectors())

print("Creating system...")
system = ff.createSystem(
    modeller.topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1 * unit.nanometer,
    constraints=app.HBonds,
    rigidWater=True,
)

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
app.PDBFile.writeFile(simulation.topology, positions, open(f"{replicate}/minimized.pdb", "w"))
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
    app.CheckpointReporter(f"{replicate}/checkpoint.chk", int(10 * unit.nanoseconds / timestep))
)

simulation.step(int(runtime / timestep))
