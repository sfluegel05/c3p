"""
Classifies: CHEBI:65061 2,5-diketopiperazines
"""
"""
Classifies: 2,5-diketopiperazines (Any piperazinone that has a piperazine-2,5-dione skeleton.)

This function determines if the input molecule contains a 2,5-diketopiperazine core.
The core is defined as a six-membered ring with the following sequence of atoms:
    N – C(=O) – C – N – C(=O) – C
where the carbonyl carbons (CX3(=O)) are at positions 2 and 5 of the ring.
The SMARTS pattern is designed to require that the atoms belong to a ring of size 6.
This helps to avoid false positives from related but non‐diketopiperazine scaffolds (e.g.
5‑membered imidazolidine-2,4-diones).
"""

from rdkit import Chem

def is_2_5_diketopiperazines(smiles: str):
    """
    Determines if a molecule is a 2,5-diketopiperazine based on its SMILES string.
    A 2,5-diketopiperazine is defined as any piperazinone that has a piperazine-2,5-dione
    skeleton. The core structure consists of a six-membered ring with two nitrogen atoms at
    positions 1 and 4, and two carbonyl groups at positions 2 and 5 (i.e.
    N - C(=O) - C - N - C(=O) - C).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2,5-diketopiperazine, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the 2,5-diketopiperazine core.
    # We require a six-membered ring with this connectivity:
    # [N;R6] - [CX3](=O) - [CX4;R6] - [N;R6] - [CX3](=O) - [CX4;R6]
    # This pattern excludes similar cores found in 5-membered rings.
    diketopiperazine_smarts = "[N;R6][CX3](=O)[CX4;R6][N;R6][CX3](=O)[CX4;R6]"
    pattern = Chem.MolFromSmarts(diketopiperazine_smarts)
    if pattern is None:
        return False, "Error in SMARTS pattern definition"
    
    # Check if the molecule contains the diketopiperazine substructure
    if mol.HasSubstructMatch(pattern):
        return True, "Molecule contains a piperazine-2,5-dione (2,5-diketopiperazine) skeleton"
    else:
        return False, "Molecule does not contain the required piperazine-2,5-dione skeleton"

# Example usage:
# Uncomment the following lines to test:
# test_smiles = [
#     "O=C1NC(=O)CN1",  # simplest case, should be True
#     "CC1C(=O)NC(=O)N1",  # 5-methylimidazolidine-2,4-dione, 5-membered ring, should be False
# ]
# for smi in test_smiles:
#     result, reason = is_2_5_diketopiperazines(smi)
#     print(smi, result, reason)