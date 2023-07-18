import os
import re

# Get the directory path of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Find all .out files in the script's directory
out_files = [filename for filename in os.listdir(script_dir) if filename.endswith('.out')]

# Reference values for chemical shifts
reference_values = {
    'C': 195.4494,
    'H': 32.7715
}

# Process each .out file
for out_file in out_files:
    file_path = os.path.join(script_dir, out_file)

    # Open the file for reading
    with open(file_path, 'r') as file:
        data = file.read()

    # Find all matches of atom information using regular expressions
    matches = re.findall(r'Atom:\s+(\d+)\s+(\w+)\s+.*?isotropic\s+=\s+([\d.-]+)', data, re.DOTALL)

    # Iterate over the matches and calculate the chemical shift
    for match in matches:
        atom_number = match[0]
        atom_name = match[1]
        isotropic_value = float(match[2])

        if atom_name in reference_values:
            reference_value = reference_values[atom_name]
            chemical_shift = reference_value - isotropic_value
            print(f'Atom {atom_number} {atom_name}: chemical shift = {chemical_shift:.4f}')
        elif atom_name == 'H':
            chemical_shift = 32.7715 - isotropic_value
            print(f'Atom {atom_number} {atom_name}: chemical shift = {chemical_shift:.4f}')

