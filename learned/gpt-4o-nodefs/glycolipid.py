"""
Classifies: CHEBI:33563 glycolipid
"""
from rdkit import Chem

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    A glycolipid has a sugar moiety linked to a lipid component.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a glycolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more accurate pattern for sugar moieties - pyranose (6-member) or furanose (5-member) rings
    sugar_patterns = [
        Chem.MolFromSmarts("C1(CO)OC(O)C(O)C1O"),  # glucose pyranose form
        Chem.MolFromSmarts("C1OCC(O)C(O)C1O"),    # furanose potential base
    ]

    # Check for any sugar moiety
    sugar_found = any(mol.HasSubstructMatch(pattern) for pattern in sugar_patterns)
    if not sugar_found:
        return False, "No sugar moiety found"

    # Check for long carbon chains representing lipids; ensure more than just a simple 'CCCC' pattern
    lipid_pattern = Chem.MolFromSmarts("C(CCCCCCCCCCCCCCCCCCCCC)")

    if not mol.HasSubstructMatch(lipid_pattern):
        return False, "No long carbon chain (lipid) found"

    return True, "Contains sugar moiety linked to a lipid component"