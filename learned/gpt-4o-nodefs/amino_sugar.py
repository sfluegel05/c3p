"""
Classifies: CHEBI:28963 amino sugar
"""
from rdkit import Chem

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    Amino sugars are typically sugars where a hydroxyl group is replaced by an amino group,
    often with an N-acetamido modification.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an amino sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define patterns for general sugar rings and amino group modifications
    pyranose_pattern = Chem.MolFromSmarts("[C@H1]1COC(O)C(O)C(O)C1") 
    furanose_pattern = Chem.MolFromSmarts("[C@H1]1CCOC(O)C1")
    amino_replacement_pattern = Chem.MolFromSmarts("[C@H](N)[C@H]")
    acetamido_replacement_pattern = Chem.MolFromSmarts("[N][C@H]([CH3])C(=O)")

    # Check for sugar moiety (detecting pyranose or furanose structures)
    has_sugar_ring = mol.HasSubstructMatch(pyranose_pattern) or mol.HasSubstructMatch(furanose_pattern)
    if not has_sugar_ring:
        return False, "No sugar ring structure detected"

    # Check for either amino group or N-acetylation modification
    has_amino_modification = mol.HasSubstructMatch(amino_replacement_pattern) or mol.HasSubstructMatch(acetamido_replacement_pattern)
    if not has_amino_modification:
        return False, "No amino or acetamido group modification found"

    return True, "Contains sugar ring with an appropriate amino modification, consistent with an amino sugar"