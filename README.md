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

