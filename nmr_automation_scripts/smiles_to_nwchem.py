from rdkit import Chem
from rdkit.Chem import AllChem
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

# Example usage:
if len(sys.argv) < 3:
    print("Usage: python smiles_to_nwchem.py [SMILES] [output_file]")
else:
    smiles = sys.argv[1]
    output_file = sys.argv[2]
    generate_3d_coordinates(smiles, output_file)
