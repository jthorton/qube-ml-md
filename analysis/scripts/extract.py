import click
import MDAnalysis as mda


@click.command()
@click.argument("topology", type=click.Path(exists=True))
@click.argument("trajectory", type=click.Path(exists=True))
@click.argument("output_name", type=click.STRING)
@click.argument("resname", type=click.STRING)
def main(topology: str, trajectory: str, output_name: str, resname: str):
    """
    Extract the solute and a sphere of solvent with in a 27 angstrom radius from the trajectory.

    TOPOLOGY:
        The name of the pdb with the full system topology
    TRAJECTORY:
        The name of the dcd file which contains the trajectory of the full system.
    OUTPUT_NAME:
        The name that should be used for the output files with no suffix.
    RESNAME:
        The resiude name of the solute used to extract it.
    """

    # load in the pdb and dcd
    u = mda.Universe(topology, trajectory)
    # cut out a small sphere around the solute molecule
    zone = u.select_atoms(f"sphzone 27.0 ( resname {resname} )").residues
    # write out the molecule to file
    zone.atoms.write(f"{output_name}.pdb", frames=[0])
    # write out the rest of traj to file as dcd
    zone.atoms.write(f"{output_name}.dcd", frames="all")


if __name__ == "__main__":
    main()
