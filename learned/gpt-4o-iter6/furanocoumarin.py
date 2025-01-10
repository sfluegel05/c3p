"""
Classifies: CHEBI:24128 furanocoumarin
"""
from rdkit import Chem

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    A furanocoumarin is defined as a furan ring fused with a coumarin structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a furanocoumarin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define substructure patterns for furan and coumarin
    furan_pattern = Chem.MolFromSmarts('c1ccoc1')  # A simple furan ring
    coumarin_pattern = Chem.MolFromSmarts('O=c1cc2ccccc2oc1')  # A common coumarin structure

    # Check for the presence of a furan ring
    furan_matches = mol.GetSubstructMatches(furan_pattern)
    if not furan_matches:
        return False, "No furan-like ring found"

    # Check for the presence of a coumarin structure
    coumarin_matches = mol.GetSubstructMatches(coumarin_pattern)
    if not coumarin_matches:
        return False, "No coumarin-like structure found"

    # Verify fusion by checking atom connectivity between matched furan and coumarin substructures
    for furan_match in furan_matches:
        for coumarin_match in coumarin_matches:
            # Check if there is a direct atom connection between end atoms of furan and coumarin rings
            if any(mol.GetBondBetweenAtoms(furan_atom, coumarin_atom) for furan_atom in furan_match for coumarin_atom in coumarin_match):
                return True, "Furan ring and coumarin structure are fused"

    return False, "Furan and coumarin structures not properly fused"