from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess
import time
import sys

def generate_3d_coordinates(smiles, output_file):
    molecule = Chem.MolFromSmiles(smiles)
    molecule = Chem.AddHs(molecule)
    AllChem.EmbedMolecule(molecule)
    AllChem.UFFOptimizeMolecule(molecule)

    with open(output_file, 'w') as f:
        f.write("start\n\n")
        f.write('title "{}"\n\n'.format(smiles))
        f.write("echo\n")
        f.write("geometry {}\n".format(smiles))

        for atom in molecule.GetAtoms():
            position = molecule.GetConformer().GetAtomPosition(atom.GetIdx())
            x, y, z = position.x, position.y, position.z
            line = "   {:<2s} {:>10.4f} {:>10.4f} {:>10.4f}\n".format(atom.GetSymbol(), x, y, z)
            f.write(line)
        f.write("end\n\n")
        f.write("\nbasis spherical\n")
        f.write("* library 6-311G\n")
        f.write("end\n\n")
        f.write("dft\n")
        f.write("direct\n")
        f.write("xc pbe0\n")
        f.write('noprint "final vectors analysis" multipole\n')
        f.write("end\n\n")
        f.write("set geometry {}\n".format(smiles))
        f.write("task dft optimize\n\n")
        f.write("property\n")
        f.write("   shielding\n")
        f.write("end\n\n")
        f.write("task dft property\n")

def run_nwchem(commands):
    print("Running NWChem...")
    start_time = time.time()

    subprocess.run(commands, shell=True)

    end_time = time.time()
    execution_time = end_time - start_time
    print("NWChem execution time: {:.2f} seconds".format(execution_time), file=sys.stderr)

# Check if SMILES input is provided as a command-line argument
if len(sys.argv) > 1:
    # Use the command-line argument as the SMILES input
    smiles = sys.argv[1]
else:
    # Prompt for SMILES input
    smiles = input("Enter SMILES: ")

# Generate NWChem input and output file names
output_file = "{}.nw".format(smiles)
output_file_name = "{}.out".format(smiles)

generate_3d_coordinates(smiles, output_file)

# Build NWChem command
nwchem_command = "nwchem {} > {}".format(output_file, output_file_name)

# Build final command with additional options
final_command = "mpiexec --use-hwthread-cpus -np 1 {}".format(nwchem_command)

# Run NWChem with the provided command
run_nwchem(final_command)

