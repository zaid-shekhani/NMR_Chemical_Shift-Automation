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
end

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

## Explanation:

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

### A. Preamble


```none
start

title "ethanol"

echo
```

The `start` statement indicates the beginning of the input file's execution. It is used to initialize the NWChem job and set the stage for defining the various sections that follow.

The `title` statement is used to provide a descriptive title or label for the NWChem job. In this case, the title is set to "ethanol," indicating that the job involves calculations related to the molecule ethanol.

The `echo` section allows users to review and verify the input parameters used in the calculation, which can be helpful for debugging or documentation purposes.

### B. Molecular Geometry

```none
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
end
```

This section defines the molecular geometry of the system under study. It specifies the atomic symbols (Oxygen, Carbon, Hydrogen) along with their respective Cartesian coordinates. The molecule described here is ethanol, with its atoms placed in 3D space.

### C. Basis Set


```none
basis spherical
* library 6-311G
end
```

This section sets the basis set used for the electronic structure calculations. A basis set is a set of mathematical functions used to describe how electrons are distributed around atoms in molecules. It helps approximate the complex behavior of electrons during computational chemistry calculations. The `6-311G` basis set is specified here, and the `spherical` keyword indicates that the basis functions are spherically symmetric.

### D. Density Functional Theory (DFT)

```none
dft
direct
xc pbe0
noprint "final vectors analysis" multipole
end
```

In this section, the Density Functional Theory (DFT) settings for the calculation are specified. DFT is a quantum mechanical method used to study electronic properties of molecules. It is used to study the electronic structure and properties of molecules and materials. It focuses on the electron density, representing the probability of finding electrons at specific locations around atomic nuclei. By analyzing electron distribution, DFT enables us to predict molecular properties, optimize geometries, and investigate chemical processes with a good balance between accuracy and computational efficiency. The keywords used here are as follows:

- `direct`: This keyword tells NWChem to use a direct SCF (Self-Consistent Field) algorithm to optimize the molecular geometry and calculate properties.
- `xc pbe0`: The exchange-correlation functional to be used is specified as pbe0. PBE0 is a hybrid functional that combines the Perdew-Burke-Ernzerhof (PBE) exchange functional with a portion of Hartree-Fock exchange.
- `noprint "final vectors analysis" multipole`: These keywords suppress printing of the final vectors analysis and multipole expansion results during the calculation.

### E. Geometry Optimization Task

```none
set geometry ethanol
task dft optimize
```

This section sets the geometry for the current calculation to the previously defined ethanol geometry. The task specified here is to perform a geometry optimization using DFT (`task dft optimize`). The optimization will find the molecular geometry that corresponds to the minimum energy on the potential energy surface.

### F. Property Calculation

```none
property
  shielding
end
task dft property
```

In this section, a property calculation is performed. The keyword `property` specifies that the properties of the system are to be calculated. The specific property calculated here is the nuclear magnetic shielding tensor (shielding) for each nucleus in the molecule. The shielding tensor provides information about the local electronic environment around each nucleus. The `task dft property` command executes the property calculation.

In summary, this NWChem input file defines the molecular geometry of ethanol, sets up the basis set and DFT settings, performs a geometry optimization, and calculates the nuclear magnetic shielding tensor as a property of the optimized ethanol molecule using the PBE0 exchange-correlation functional.
User



## Generating NWChem Input File from SMILES

To make things easier, I have created a Python script that generates the NWChem input file by providing the SMILES representation of the molecule. You can find the Python Script (`smiles_to_nwchem.py`) in this repository's [nmr_automation_scripts](./nmr_automation_scripts) folder.

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

2. Save the Python script `smiles_to_nwchem.py`, which you can find in this repository's [nmr_automation_scripts](./nmr_automation_scripts) folder, in the 'CCO' folder.

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

3. The calculation may take some time depending on the complexity of the system and the computational resources available. In my case, with an Intel® Core™ i7-8565U CPU @ 1.80GHz × 8 processor and 16 GB of RAM, it took approximately 2 minutes and 56 seconds. Please remain patient while the calculation runs.

4. Once the calculation is complete, check the main directory ("CCO"), which should now contain multiple files. Our main interest is the `ethanol.out` file, which contains the output of the NWChem calculation.

**Note:** I have provided a sample `ethanol.out` file in the Resources folder of this repository, along with the `ethanol.nw` input file. You can refer to these files for reference and comparison.

By following these steps, you can run NWChem and perform calculations using the provided input file. The output file (`ethanol.out`) will contain the results and information generated by NWChem for further analysis and interpretation.

### Analysis of Output and Calculation of NMR Chemical Shifts

To calculate the <sup>1</sup>H and <sup>13</sup>C NMR chemical shift values for ethanol, we need the output file for TMS as a reference. If you haven't done so already, follow the steps described earlier to generate the NWChem output file (`tms.out`) for TMS. If needed, you can find both the input and output files in the Resources folder of this repository.

Now let's discuss some important points about ethanol's chemical environment. Ethanol's NMR environment refers to the local electronic and magnetic environment surrounding its hydrogen nuclei. It has three different hydrogen environments: the methyl group (-CH3), the methine group (-CH2-), and the hydroxyl group (-OH). The methyl group gives a triplet signal, the methine group gives a quartet signal, and the hydroxyl group gives a broad singlet due to hydrogen bonding. Factors such as solvent, temperature, and pH can affect the observed NMR signals.

![Image1](Images/graph.svg)
![Ethanol NMR Spectra](Images/graph(1).svg)
![Ethanol NMR Spectra](Images/graph(2).svg)

In order to calculate the chemical shift values, we need the individual **isotropic** values of each atom from the output file. Open the `tms.out` file and search for the "**Chemical Shielding Tensors (GIAO, in ppm)**" section. In my case, this section starts at line 3394. Scroll down a bit, and at line 3468, you'll find "**Atom:    2  C**" which is of interest to us. The calculated isotropic value for this carbon is **195.4494** (check line 3484). Similarly, for the hydrogen of TMS, the isotropic value we obtain is **32.7715**. You will notice that all the other hydrogens and carbons have similar isotropic values. This is because all the carbons and hydrogens in TMS are in the same chemical environment, making it an ideal reference compound.

From the isotropic values (in ppm) of TMS, we obtain:
   ```
   C: 195.4494
   H: 32.7715
   ```

Similarly, for ethanol, the isotropic values (in ppm) we obtain are:
   ```
   C: 132.5255, 174.7759
   H: 32.8701, 28.9842, 28.8955, 31.8083, 31.5748, 31.2364
   ```

Now we can calculate the chemical shift values for the two carbons of ethanol as follows:

Chemical shift = Isotropic value of carbon in TMS - Isotropic value of carbon in ethanol

Thus, the chemical shift for the first carbon of ethanol is:
   ```
   Chemical shift = 195.4494 - 132.5255 = 62.9239
   ```
And for the second carbon:You can find the chemical shift results saved in the result.txt file in your working directory.
   ```
   Chemical shift = 195.4494 - 174.7759 = 20.6735
   ```

These values are in good agreement with the experimental chemical shift values of 57.4 and 16.4, respectively.

Similarly, for the three methyl hydrogens of ethanol, we obtain the following chemical shift values: 0.9632, 1.1967, and 1.5351. Taking the average, we get 1.231666667, which is very close to the experimental value of 1.187.

For the two methine hydrogens of ethanol, we obtain the following chemical shift values: 3.7873 and 3.876. Averaging them, we get 3.83165, which is very close to the experimental value of 3.314.

**<sup>13</sup>C NMR Chemcial Shifts for Ethanol:**  
| Calculated Value    | Experimental Value | 
|----------|-----|
| 62.9239    | 57.4  |
| 20.6735   | 16.4  |

**<sup>1</sup>H NMR Chemcial Shifts for Ethanol:**  
| Calculated Value    | Experimental Value | 
|----------|-----|
| 1.231666667    | 1.187  |
| 3.83165   | 3.314  |

By calculating the chemical shift values, we can compare them to experimental values and gain insights into the chemical environment of the molecule. These values demonstrate the accuracy and usefulness of NMR calculations in predicting and interpreting NMR spectra for organic molecules.

## Automation

I have developed few Python and shell scripts(which you can find in this repository's [nmr_automation_scripts](./nmr_automation_scripts) folder) to automate the manual processes involved in calculating NMR chemical shift values. With these automation scripts, you can now obtain chemical shifts with just a line of code providing the SMILES representation of the molecule. This automation significantly reduces the manual work required compared to the ISiCLE version.

### Scripts Overview

1. **auto1.py**: This Python script converts SMILES notation to NWChem input files. It takes the SMILES representation of the molecule as input and generates the corresponding NWChem input file (with a `.nw` extension).

2. **auto2.py**: This Python script extracts and calculates chemical shifts from NWChem output files. It reads the NWChem output file and calculates the chemical shift values for each atom type in the molecule.

3. **run_scripts.sh**: This shell script integrates the two Python scripts, `auto1.py` and `auto2.py`, into a seamless automation process. By running this shell script with the appropriate arguments, you can automate the entire calculation of NMR chemical shifts.

### Automation Steps

Before running the automation, ensure you have completed the prerequisites mentioned in the [Prerequisites](#prerequisites) section.

1. Activate the `rdkit-env` environment by running the following command in your terminal:
   ```
   conda activate rdkit-env
   ```
2. Make the shell script executable by running the following command in your terminal:
   ```
   chmod +x run_scripts.sh
   ```
3. Now, you can run the shell script with the desired SMILES representation to obtain NMR chemical shifts. For example, to calculate chemical shifts for the molecule with SMILES "CCO," execute the following command:
   ```
   ./run_scripts.sh "CCO"
   ```

By following these steps, the automation will convert the SMILES to the NWChem input file and perform the necessary calculations to obtain the chemical shift values. The results will be displayed in the terminal.

![Image2](Images/Screenshot.png)

You can find the chemical shift results saved in the **result.txt** file in your working directory.

**Note:** Make sure to replace "CCO" with the appropriate SMILES representation for the molecule you want to calculate the chemical shifts for.

Feel free to use these automation scripts to expedite your NMR chemical shift calculations and streamline the process.
