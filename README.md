# NWChem Tutorial: Calculating NMR Chemical Shift Values

## Overview
Welcome to the NWChem Tutorial for calculating Nuclear Magnetic Resonance (NMR) chemical shift values! In this tutorial, I will explore how NWChem, a powerful computational chemistry software suite, can be utilized to calculate NMR chemical shifts of chemical compounds.

## Introduction
Nuclear Magnetic Resonance (NMR) spectroscopy is a valuable technique for chemical identification in complex unknown samples. Traditionally, researchers rely on databases containing chemical shift spectra to match experimental NMR data with known compounds. However, building comprehensive libraries of authentic standards experimentally can be challenging and time-consuming.

To overcome these limitations, computational chemistry approaches offer a practical solution. By utilizing NWChem, a widely-used software suite for quantum chemistry calculations, we can perform in silico calculations to determine molecular attributes such as NMR chemical shifts. This enables us to expand reference libraries and enhance compound identification.

## Installation
The following tutorial provides instructions for installing NWChem on Linux. If you are using a different operating system, such as Windows or macOS, please visit the [NWChem website](https://nwchemgit.github.io/Download.html#nwchem-availability-in-linux-distributions) for specific instructions tailored to your OS.
### Installation Steps
1. Open the terminal.
2. Update your package list by running the following command:
   ```shell
   sudo apt update
3. Install NWChem by running the following command:
   ```shell
   sudo apt install nwchem
4. Verify the installation by running the following command to check the NWChem version:
   ```shell
   nwchem --version
This should display the version information if NWChem is installed correctly.

5. Congratulations! You have successfully installed NWChem on your system.
## Input Files

NWChem uses input files to specify the details of the calculations or simulations you want to perform. The input file is a text file that contains various sections and keywords to define the molecular system, computational settings, and desired calculations. Let's take a look at an example input file:

```none
start

title "ethanol"

echo
geometry ethanol
   O       1.4277    -0.5650    -0.3289
   C       0.4703     0.4443    -0.1587
   C      -0.9137    -0.1644     0.0308
   H       1.6014    -0.9462     0.5713
   H       0.4633     1.0865    -1.0645
   H       0.7274     1.0839     0.7160
   H      -0.9322    -0.7928     0.9464
   H      -1.6672     0.6445     0.1353
   H      -1.1770    -0.7909    -0.8477

basis spherical
* library 6-311G
end

dft
direct
xc pbe0
noprint "final vectors analysis" multipole
end

set geometry ethanol
task dft optimize
property
   shielding
end
task dft property
```

**Explanation:**

- The `start` keyword marks the beginning of the input file.
- The `title` keyword specifies a title for the calculation.
- The `echo` keyword enables the printing of input sections and commands to the output file.
- The `geometry` section defines the molecular geometry of the ethanol molecule in Cartesian coordinates.
- The `basis` section specifies the basis set to be used, which is the 6-311G basis set in this case.
- The `dft` section sets up the DFT calculation using the PBE0 exchange-correlation functional and enables direct calculation methods.
- The `noprint` keyword suppresses printing of the final vectors analysis and multipole moments.
- The `set` keyword sets the geometry to the previously defined ethanol geometry.
- The `task dft optimize` line instructs NWChem to perform a DFT optimization of the molecular structure.
- The `property` section specifies the property calculation to be performed, which is the calculation of NMR shielding tensors in this case.
- The final `task dft property` line instructs NWChem to perform the property calculation.

## Generating NWChem Input File from SMILES

To make things easier, I have created a Python script that generates the NWChem input file by providing the SMILES representation of the molecule. You can find the Python Script (`smiles_to_nwchem.py`) in this repository's `nmr_automation_scripts` folder.

### Prerequisites

Before using the Python script, you need to ensure that RDKit is installed. RDKit is a powerful open-source cheminformatics library that provides various tools for molecular analysis and manipulation. Here's a step-by-step guide to installing RDKit using Anaconda:

1. Install Anaconda by following the instructions provided on the [Anaconda website](https://www.anaconda.com/products/individual).

2. Once Anaconda is installed, open a terminal and create a new environment named `rdkit-env` by running the following command:
   ```
   conda create -n rdkit-env -c conda-forge rdkit
   ```
3. Activate the `rdkit-env` environment using the following command:
   ```
   conda activate rdkit-env
   ```

This will ensure that you are working within the environment where RDKit is installed.

### Generating NWChem Input File

To use the Python script and convert the SMILES representation to an NWChem input file, follow these steps:

1. Create a new folder, let's name it 'CCO', where you will save the generated input files.

2. Save the Python script `smiles_to_nwchem.py`, which you can find in this repository's nmr_automation_scripts folder, in the 'CCO' folder.

3. Open the terminal and navigate to the 'CCO' folder directory.

4. Activate the `rdkit-env` environment by running the following command:
   ```
   conda activate rdkit-env
   ```

This step ensures that the script has access to the RDKit library.

5. Run the Python script with the desired SMILES representation to generate the NWChem input file. For example, to create an input file for ethanol, execute the following command:
   ```
   python smiles_to_nwchem.py "CCO" ethanol.nw
   ```

In this case, the SMILES representation "CCO" is provided as the first argument, and the output NWChem input file is named `ethanol.nw`. You can replace "CCO" with the appropriate SMILES for the molecule you want to generate an input file for.

**Note:** Make sure to always save the NWChem input files with the `.nw` extension at the end.

By following these steps, you can use the Python script to automatically generate NWChem input files based on the SMILES representation of the molecule you specify. This simplifies the process and saves you time when setting up calculations for different molecules.

### Running NWChem

Now that we have our NWChem input file (`ethanol.nw`), let's proceed to run NWChem and perform the calculations. Follow these steps:

1. Open a terminal and navigate to the main directory, "CCO," if you are not already there.

2. Run the following command to execute NWChem with the input file:
   ```
   time nwchem ethanol.nw > ethanol.out 2> ethanol.log
   ```

This command runs NWChem using the `ethanol.nw` input file, directs the standard output to the `ethanol.out` file, and redirects the error output to the `ethanol.log` file.

**Note:** The `time` command is used to measure the execution time of the NWChem calculation. It is optional but provides information on the runtime.

3. The calculation may take some time depending on the complexity of the system and the computational resources available. In my case, with an Intel® Core™ i7-8565U CPU @ 1.80GHz × 8 processor and 16 GB of RAM, it took approximately 3 minutes and 9 seconds. Please remain patient while the calculation runs.

4. Once the calculation is complete, check the main directory ("CCO"), which should now contain multiple files. Our main interest is the `ethanol.out` file, which contains the output of the NWChem calculation.

**Note:** I have provided a sample `ethanol.out` file in the Resources folder of this repository, along with the `ethanol.nw` input file. You can refer to these files for reference and comparison.

By following these steps, you can run NWChem and perform calculations using the provided input file. The output file (`ethanol.out`) will contain the results and information generated by NWChem for further analysis and interpretation.


