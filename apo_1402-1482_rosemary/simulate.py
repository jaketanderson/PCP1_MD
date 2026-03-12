import os
import sys

import openmm
from openff import toolkit
from openff.units import unit as openff_unit
from openmm import app, unit

from openff.pablo import STD_CCD_CACHE, ResidueDefinition, topology_from_pdb
from openff.toolkit import ForceField
from openff.interchange.interop.openmm import to_openmm_positions

replicate = int(sys.argv[1])
print(f"This is replicate #{replicate}.")

timestep = 2 * unit.femtosecond
runtime = 1050 * unit.nanoseconds
temperature = 298.15 * unit.kelvin

topology = topology_from_pdb("PCP1_solvated.pdb")

ff = ForceField(
    "openff_no_water-3.0.0-alpha0.offxml",
    "opc-1.0.2.offxml",
)

interchange = ff.create_interchange(topology)

print("Creating system...")

interchange["vdW"].cutoff = 1.0 * openff_unit.nanometer
interchange["Electrostatics"].cutoff = 1.0 * openff_unit.nanometer
system = interchange.to_openmm_system()

with open(f"workspace/{replicate}/system.xml", "w") as f:
    f.write(openmm.XmlSerializer.serialize(system))

integrator = openmm.LangevinMiddleIntegrator(temperature, 1 / unit.picosecond, timestep)
# Use the int `replicate` as the random number seed
integrator.setRandomNumberSeed(replicate)
platform = openmm.Platform.getPlatformByName("CUDA")
print("CUDA_VISIBLE_DEVICES:", os.environ.get("CUDA_VISIBLE_DEVICES"))
properties = {"Precision": "mixed"}

simulation = app.Simulation(
    interchange.to_openmm_topology(),
    system,
    integrator,
    platform,
    properties,
)

simulation.context.setPositions(to_openmm_positions(interchange))

print("Minimizing energy...")
simulation.minimizeEnergy(tolerance=2.5 * unit.kilojoules_per_mole / unit.nanometer)
positions = simulation.context.getState(positions=True).getPositions()
with open(f"workspace/{replicate}/minimized.pdb", "w") as f:
    app.PDBFile.writeFile(simulation.topology, positions, f)
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
