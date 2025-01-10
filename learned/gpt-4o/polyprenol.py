"""
Classifies: CHEBI:26199 polyprenol
"""
from rdkit import Chem

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    A polyprenol has more than one isoprene unit and ends with an alcohol group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Broader pattern for isoprene units - consider variations, ensures two C=C bonds
    isoprene_pattern_a = Chem.MolFromSmarts("C=C-C")
    isoprene_pattern_b = Chem.MolFromSmarts("C=C-C=C")
    isoprene_matches_a = mol.GetSubstructMatches(isoprene_pattern_a)
    isoprene_matches_b = mol.GetSubstructMatches(isoprene_pattern_b)
    isoprene_matches = len(isoprene_matches_a) + len(isoprene_matches_b)
    
    # Check if there are more than one isoprene units across both patterns
    if isoprene_matches < 2:
        return False, f"Found {isoprene_matches} isoprene units, need more than one"

    # Ensure terminal alcohol group - check for groups bonded to terminal carbons
    alcohol_pattern = Chem.MolFromSmarts("[CX4;H2,O]O")
    terminal_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetDegree() == 1]
    if not any(mol.GetSubstructMatches(alcohol_pattern, atomId=idx) for idx in terminal_atoms):
        return False, "No terminal alcohol group found or alcohol group not at terminal position"

    return True, "Contains more than one isoprene unit and a terminal alcohol group"