# qube-ml-md
A collection of scripts to run and analyse MD simulations powered by QUBE/ML force fields.

## Installation

### 1. Clone the repository:

```
git clone https://github.com/jthorton/qube-ml-md.git
cd qube-ml-md
```

### 2. Install and activate the conda environment:

```
conda env create -f env.yml
conda activate qube-ml-md
```

### 3. MD scripts

First copy the `run_md.py` script into the working directory which should contain a pdb topology file and solute and solvent 
force field files in `OpenMM` xml format. 

#### MM simulations

As an example we will run a pure MM seed/production simulation using the parsley force fields and coumarin153 topology provided.

To see the required inputs and options of the script run
```
python run_md.py --help
```

We see that the script at minimum requires an input topology file force fields and the number of seed snapshots that should be saved.
The following command will run the pure MD simulations using the CUDA platform which is available with an NVIDIA GPU.

```
python run_md.py mixture.pdb coumarin_parsley.xml acetonitrile_parsley.xml 2 -p CUDA
```

This will spawn a seed simulation which will save 2 snapshots seperated by 50ps of MD. Each seed state is then used for a production simulation.
The production simulations will be run in separate folders created where the script is run.

#### MM/ML Simulations

We can also use the script to create a hybrid system using ani2x. To do this we have to add a flag to tell the script we want a hybrid 
system and then tell it the residue name of the solute molecule we want to apply ani2x to. Note the platforms `CUDA` or `OpenCL` must be used here.

```
python run_md.py mixture.pdb coumarin_parsley.xml acetonitrile_parsley.xml 2 -p CUDA -ani -solute COU
```

### 4. Extract

Once the production simulation is complete a script is provided in the `analysis` folder to extract a trajectory for the solute 
enclosed by a sphere of solvent within 27 angstroms of the centre of the solute. First move the script to the production folder 
and copy in the starting topology. Then run 

```
python extract.py mixture.pdb trajectory.dcd sphere COU
```

for help with this script try:

```
python extract.py --help
```
