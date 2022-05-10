from simtk.openmm import app
import click


@click.command()
@click.argument("topology_file", type=click.Path(exists=True))
@click.argument("solute_xml", type=click.Path(exists=True))
@click.argument("solvent_xml", type=click.Path(exists=True))
@click.argument("output_file", type=click.STRING)
def add_sites(topology_file: str, solute_xml: str, solvent_xml:str, output_file: str):
    """
    Add vsites described in the solute xml file to the topology file and write out the new positions to the output file.
    
    TOPOLOGY_FILE:
        The name of the pdb file that contains the topology of the system.

    SOLUTE_XML:
        The force field xml file used to parametersise the solute which contains the vsite definitions.

    SOLVENT_XML:
        The force field xml file used to parameterise the solvent which contains any vsite definitions.

    OUTPUT_FILE:
        The name of the file which the final topology should be written to including vsites.
    """

    pdbfile = app.PDBFile(topology_file)
    ff = app.ForceField(solute_xml, solvent_xml)
    modeller = app.Modeller(topology=pdbfile.topology, positions=pdbfile.positions)
    modeller.addExtraParticles(forcefield=ff)
    app.PDBFile.writeFile(topology=modeller.topology, positions=modeller.positions, file=open(output_file, "w"))


if __name__ == "__main__":
    add_sites()
