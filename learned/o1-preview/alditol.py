"""
Classifies: CHEBI:17522 alditol
"""
"""
Classifies: CHEBI:18385 alditol
"""

from rdkit import Chem

def is_alditol(smiles: str):
    """
    Determines if a molecule is an alditol based on its SMILES string.
    An alditol is an acyclic polyol with the general formula HOCH2[CH(OH)]nCH2OH,
    formally derivable from an aldose by reduction of the carbonyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alditol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Remove the check for cyclic molecules to allow alditol chains within larger structures
    # Define the maximum length of the alditol chain to search for
    max_n = 10  # Adjust as needed

    # Initialize a flag to indicate if an alditol chain is found
    alditol_found = False

    # Try matching alditol chains of varying lengths
    for n in range(1, max_n + 1):  # n is the number of internal CH(OH) units
        # Build the SMARTS pattern for the alditol chain
        chain = '-'.join(['[C@H](O)'] * n)
        pattern_str = f'[OH]-[CH2]-{chain}-[CH2]-[OH]'
        pattern = Chem.MolFromSmarts(pattern_str)
        if pattern is None:
            continue  # Skip invalid patterns

        # Search for the pattern in the molecule
        matches = mol.GetSubstructMatches(pattern)
        if matches:
            alditol_found = True
            return True, f"Molecule contains an alditol chain of length {n + 2} carbons"

    # If no alditol chain is found
    if not alditol_found:
        return False, "Molecule does not contain an alditol chain"