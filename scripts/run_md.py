from simtk.openmm import app
from simtk import openmm, unit
from typing import Optional, Tuple
from openmmml import MLPotential
from typing_extensions import Literal
import click
import os


def build_simulation(
    system: openmm.System,
    topology: app.Topology,
    temperature: unit.Quantity,
    platform: str,
) -> app.Simulation:
    """
    Build the openMM simulation from the openmm system.

    Args:
        system:
            The openmm system that should be used to build the simulation object.
        topology:
            The openmm topology which represents the system.
        temperature:
            The temperature the simulation should be run at.
        platform:
            The target platform the simulation should be run on.
    """
    time_step = 0.5 * unit.femtoseconds  # simulation timestep
    friction = 1 / unit.picosecond  # collision rate
    integrator = openmm.LangevinIntegrator(temperature, friction, time_step)
    # get the CUDA platform
    platform = openmm.Platform.getPlatformByName(platform)
    if platform == "CUDA" or platform == "OpenCL":
        properties = {"Precision": "double"}
    else:
        properties = {}
    # Set up an OpenMM simulation.
    simulation = openmm.app.Simulation(
        topology=topology,
        system=system,
        integrator=integrator,
        platform=platform,
        platformProperties=properties,
    )

    return simulation


def build_system(
    topology: app.Topology,
    positions: unit.Quantity,
    force_field: app.ForceField,
    use_ani2x: bool,
    solute_resname: Optional[str] = None,
) -> Tuple[openmm.System, app.Topology, unit.Quantity]:
    """
    Build an openmm mm system from the topology and forcefield files.

    Args:
        topology:
            The topology of the system
        positions:
            The initial poistions of the system, used to calculate the positions of any virtual sites
        force_field:
            The force field that should be used to parameterise the system.
        use_ani2x:
            If we should use ani2x for the internal energy of the solute?
        solute_resname:
            The residue name of the solute in the topology, this must be unique.
    """
    print("Building system")
    try:
        mm_system = force_field.createSystem(
            topology=topology,
            nonbondedMethod=app.PME,
            nonbondedCutoff=1 * unit.nanometer,
            ewaldErrorTolerance=1e-4,
        )
        initial_positions = positions
        initial_topology = topology

    except ValueError:
        modeller = app.Modeller(topology=topology, positions=positions)
        modeller.addExtraParticles(forcefield=force_field)
        mm_system = force_field.createSystem(
            topology=modeller.topology,
            nonbondedMethod=app.PME,
            nonbondedCutoff=1 * unit.nanometer,
            ewaldErrorTolerance=1e-4,
        )
        initial_positions = modeller.positions
        initial_topology = modeller.topology

    if use_ani2x:
        if solute_resname is None:
            raise RuntimeError(
                "The solute residue name must be supplied and be unique to use ani2x."
            )
        print("Creating a mixed MM/ML system using ani2x.")
        potential = MLPotential("ani2x")
        ml_atoms = [
            atom.index
            for atom in topology.atoms()
            if atom.residue.name == solute_resname
        ]
        mm_system = potential.createMixedSystem(topology, mm_system, ml_atoms)

    return mm_system, initial_topology, initial_positions


platforms = Literal["CPU"]


@click.command()
@click.argument("topology_file", type=click.Path(exists=True))
@click.argument("solute_xml", type=click.Path(exists=True))
@click.argument("solvent_xml", type=click.Path(exists=True))
@click.argument("seeds", type=int)
@click.option(
    "-ani",
    "--use_ani2x",
    is_flag=True,
    help="If ani2x should be used for the internal energy of the solute or not.",
)
@click.option(
    "-solute",
    "--solute_resname",
    type=click.STRING,
    help="The unique residue name of the solute needed for ani2x.",
)
@click.option(
    "-p",
    "--platform",
    default="CUDA",
    type=click.Choice(["CUDA", "OpenCL", "CPU", "Reference"]),
)
def run(
    topology_file: str,
    solute_xml: str,
    solvent_xml: str,
    seeds: int,
    use_ani2x: bool = False,
    solute_resname: Optional[str] = None,
    platform: str = "CUDA",
):
    """
    The main function of the script which loads the system and runs the seed simulations followed by the production runs.

    TOPOLOGY_FILE:
        The name of the pdb file that contains the topology of the system.

    SOLUTE_XML:
        The xml file used to parameterise the solute.

    SOLVENT_XML:
        The xml file used to parametersie the solvent.

    SEEDS:
        The number of seed snapshots to save and the number of production runs.
    """
    # load the topology and xml force field files
    pdbfile = app.PDBFile(topology_file)
    # get input positions
    positions = pdbfile.getPositions()
    ff = app.ForceField(solute_xml, solvent_xml)
    # create the system and get the correct initial positions
    system, initial_topology, initial_positions = build_system(
        topology=pdbfile.topology,
        positions=positions,
        force_field=ff,
        use_ani2x=use_ani2x,
        solute_resname=solute_resname,
    )
    # add the barostat
    temperature = 300 * unit.kelvin
    system.addForce(openmm.MonteCarloBarostat(1 * unit.atmosphere, temperature))
    print("Writing system to system.xml")
    with open("system.xml", "w") as output:
        output.write(openmm.XmlSerializer.serializeSystem(system))
    simulation = build_simulation(
        system=system,
        topology=initial_topology,
        temperature=temperature,
        platform=platform,
    )
    # set the initial positions
    simulation.context.setPositions(initial_positions)
    print("Minimising system energy...")
    simulation.minimizeEnergy()
    if not use_ani2x:
        # we can only do this if not using ani2x as we get a GIL lock bug in openmm
        # Randomize the velocities from a Boltzmann distribution at a given temperature.
        simulation.context.setVelocitiesToTemperature(temperature)

    print(f"Starting seed simulation saving {seeds} snapshots")

    # Length of the simulation 50ps.
    num_steps = 100_000  # number of integration steps to run

    # Logging options.
    trj_freq = 10_000  # number of steps per written trajectory frame
    data_freq = 1000  # number of steps per written simulation statistics

    # Configure the information in the output files.
    dcd_reporter = openmm.app.DCDReporter("trajectory.dcd", trj_freq)
    state_data_reporter = openmm.app.StateDataReporter(
        "data.csv",
        data_freq,
        step=True,
        potentialEnergy=True,
        volume=True,
        temperature=True,
        speed=True,
    )
    simulation.reporters.append(dcd_reporter)
    simulation.reporters.append(state_data_reporter)

    for i in range(seeds):
        # do 50 ps steps
        simulation.step(num_steps)
        print(f"Writing current positions for state {i}")
        simulation.saveState(f"seed_state_{i}.xml")

    print("Seed simulation complete.")

    print("Stating production simulations...")

    # define production settings
    # 16 ps
    prod_num_steps = 32_000
    # logging details
    prod_traj_freq = 4
    prod_data_freq = 1000
    home = os.getcwd()
    for i in range(seeds):
        # clear out the reporters
        simulation.reporters = []
        # load the sate to the simulation
        print(f"Loading state {i}")
        state_file = f"seed_state_{i}.xml"
        simulation.loadState(state_file)
        # move to the new production sim
        prod_dir = f"production_run_{i}"
        os.makedirs(prod_dir, exist_ok=True)
        os.chdir(prod_dir)
        # run the simulation
        dcd_reporter = app.DCDReporter("trajectory.dcd", prod_traj_freq)
        state_data = app.StateDataReporter(
            "data.csv",
            prod_data_freq,
            step=True,
            potentialEnergy=True,
            temperature=True,
            volume=True,
            speed=True,
        )
        simulation.reporters.append(state_data)
        print("Starting equ simulation ...")
        simulation.step(10_000)

        print("Starting main simulation")
        simulation.reporters.append(dcd_reporter)
        simulation.step(prod_num_steps)
        os.chdir(home)


if __name__ == "__main__":
    run()
