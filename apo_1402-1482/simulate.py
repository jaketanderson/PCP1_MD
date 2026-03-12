import sys
import os

import openmm
from openmm import app, unit

replicate = int(sys.argv[1])
print(f"This is replicate #{replicate}.")

timestep = 2 * unit.femtosecond
runtime = 1050 * unit.nanoseconds
temperature = 298.15 * unit.kelvin

pdb = app.PDBFile("5u3h_prepped.pdb")
ff = app.ForceField("amber14/protein.ff14SB.xml", "amber/opc_standard.xml")
modeller = app.Modeller(pdb.topology, pdb.positions)
modeller.addExtraParticles(ff)
modeller.addSolvent(
    ff,
    padding=1.2 * unit.nanometers,
    model="tip4pew",
    ionicStrength=0.150 * unit.molar,
    boxShape="dodecahedron",
)

print("Creating system...")
system = ff.createSystem(
    modeller.topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1 * unit.nanometer,
    constraints=app.HBonds,
)

with open(f"workspace/{replicate}/system.xml", "w") as f:
    f.write(openmm.XmlSerializer.serialize(system))

integrator = openmm.LangevinMiddleIntegrator(temperature, 1 / unit.picosecond, timestep)
# Use the int `replicate` as the random number seed
integrator.setRandomNumberSeed(replicate)
platform = openmm.Platform.getPlatformByName("CUDA")
print("CUDA_VISIBLE_DEVICES:", os.environ.get("CUDA_VISIBLE_DEVICES"))
properties = {"Precision": "mixed"}
#simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)

simulation.context.setPositions(modeller.positions)
print("Minimizing energy...")
simulation.minimizeEnergy(tolerance=2.5 * unit.kilojoules_per_mole / unit.nanometer)
positions = simulation.context.getState(positions=True).getPositions()
app.PDBFile.writeFile(
    simulation.topology, positions, open(f"workspace/{replicate}/minimized.pdb", "w")
)
simulation.saveState(f"workspace/{replicate}/minimized_state.xml")

print("Simulating production...")
simulation.reporters.append(
    app.StateDataReporter(
        f"workspace/{replicate}/production.log",
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
        f"workspace/{replicate}/production.dcd",
        reportInterval=int(20 * unit.picosecond / timestep),
        enforcePeriodicBox=True,
    )
)

simulation.reporters.append(
    app.CheckpointReporter(
        f"workspace/{replicate}/checkpoint.chk", int(10 * unit.nanoseconds / timestep)
    )
)

simulation.step(int(runtime / timestep))
