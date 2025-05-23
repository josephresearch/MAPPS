# Project README

> This is a supporting code for manuscript for "Prediction of small-molecule partitioning into biomolecular condensates from simulation" (2025) [Alina Emelianova, Pablo L. Garcia, Daniel Tan, and Jerelle A. Joseph].

> The following steps can be used to generate the coarse-grained parameters for small molecule using MAPPS algorithm from the generated all-atom (AA) and coarse-grained (CG) simulations data with GROMACS and LAMMPS. The folders contain the simulation files used to generate the data for the original publication. The code is not fully automated (yet) and assumes the output data from the PMF and slab simulations in both AA and CG resolutions is available. The demo files to generate such data are supplied, as well as the input files used to obtain the results presented in the publication.

## Table of Contents

- [System Requirements](#system-requirements)
- [Step 1: Self-Interaction Parameters for Small Molecules](#step-1-self-interaction-parameters-for-small-molecules)
    - [Folder Structure](#folder-structure)
    - [Usage](#usage)
    - [Steps](#steps)
- [Step 2: Get SM-Y Parameter: AA Partitioning and Parameter Fitting](#step-2-aa-partitioning-and-parameter-fitting)
    - [Step 2.1: AA Partitioning](#step-21-aa-partitioning)
        - [AA Partitioning Script](#aa-partitioning-script)
            - [Usage](#usage-1)
    - [Step 2.2: Fit CG to AA](#step-22-fit-cg-to-aa)
        - [Instructions](#instructions)
- [Step 3: Get Final SM-X Parameters](#step-3-get-final-sm-x-parameters)

## System Requirements

- **Operating System**: Linux or macOS
- **Software**:
    - GROMACS (version 2020 or later)
    - LAMMPS (version 3 Mar 2020 or later)
    - Python (version 3.7 or later)
    - Jupyter Notebook
- **Python Packages**:
    - numpy
    - pandas
    - matplotlib
    - scipy

## Step 1: Self-Interaction Parameters for Small Molecules

This folder is used to generate self-interaction parameters for small molecules from the pair PMF profiles.
### Folder Structure

- **compounds**: Contains subfolders named after each compound, each containing a PMF profile for self-interaction.

### Usage

The Jupyter notebook `write_params_molecules_self.ipynb` uses the PMF profiles to extract the parameters for the Wang-Frenkel potential and outputs them in a text file. These parameters should then be used in Step 2.2 when creating input files for small molecules.

Run the Jupyter notebook to get those parameters.

### Steps

1. Navigate to the `compounds` folder, create a subfolder with a molecule and save there a PMF profile obtained from GROMACS simulation.
2. Open and run the Jupyter notebook `write_params_molecules_self.ipynb`.
3. Use the generated parameters in Step 2.2 for creating input files for small molecules.

## Step 2: Get SM-Y parameter: AA Partitioning and Parameter Fitting

### Step 2.1: AA Partitioning

Navigate to the folder `Step2-fit-params_6Y/2.1.-Get_AA_partitioning`.

#### AA Partitioning Script

This script calculates the AA partitioning for a given small molecule or multiple small molecules based on the output from GROMACS simulations.

##### Usage

```bash
python get_AA_partitioning_6Y.py <DRUG>
```

Multiple small molecules can be input, separated by spaces.

### Step 2.2: Fit CG to AA

Navigate to the folder `Step2-fit-params_6Y/2.2-Fit_CG_to_AA_6Y`.

#### Instructions

1. Create the folder `YYYYYY/runs/config1_<DRUG>/mol_350/`.
2. Modify the file `YYYYYY/runs/config1_<DRUG>/mol_350/sm_ff.dat`:
   - Type in the molar mass.
   - Choose an initial guess for epsilon for `pair_coeff 9-41` (default 0.5) and type in the appropriate sigma value.
   - For `41-41`, type in parameters from self molecules PMF (from Step1-Get-SM-SM-Parameters).

3. Run the following commands:

   ```bash
   python make_runs.py
   bash YYYYYY/runs/config1_<DRUG>/submit.sh
   ```
After the simulations are finished, run the following script to obtain SM-Y parameter:

```
python find_matching_parameters.py <DRUG>
```

`<DRUG>` is the name of the small molecule you are working with. Replace `<DRUG>` with the actual small molecule name.

## Step 3: Get Final SM-X parameters

Navigate to the folder `Step3-Affinity-Score` and run:

```bash
bash create_SM-X-parameters.sh <DRUG>
```

The parameters are saved in the `params`. Copy the generated parameters and replace in the files generated in the .dat file in ./Step2-fit_params_6Y/2.2-Fit_CG_to_AA_6Y/data/CG_best_fit. These are ready-to-use parameters for a coarse-grained simulations with Mpipi force field.
