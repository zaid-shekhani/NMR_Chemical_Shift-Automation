#!/bin/bash

# Check if SMILES input is provided as a command-line argument
if [ -z "$1" ]; then
    # Prompt for SMILES input if not provided
    echo -n "Enter SMILES: "
    read -r smiles
else
    # Get the SMILES input from the command-line argument
    smiles="$1"
fi

# Save SMILES to a temporary file
echo "$smiles" > smiles.txt

# Run smiles_to_nwchem.py and redirect input from the temporary file
echo "Running NWChem..."
python smiles_to_nwchem.py < smiles.txt > /dev/null

# Extract file name from SMILES
filename=$(python -c "import re; import sys; sys.stdout.write(re.sub('[^A-Za-z0-9]+', '_', sys.argv[1]))" "$smiles")

# Run isotropic_extraction.py and capture the output
output=$(python isotropic_extraction.py "$filename.out")

# Print the captured output
echo "$output"

# Save the output to result.txt file
echo "$output" > result.txt

# Clean up the temporary file
rm smiles.txt

