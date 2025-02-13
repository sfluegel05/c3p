"""
Classifies: CHEBI:33563 glycolipid
"""
from rdkit import Chem

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    A glycolipid generally contains one or more sugar moieties linked to a lipid component.

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

    # Define broader patterns for sugar moieties
    sugar_patterns = [
        Chem.MolFromSmarts("C1(CO)OC(O)C(O)C1O"),  # Gluco-pyranose form
        Chem.MolFromSmarts("C1OC(O)C(O)C(O)C1"),  # Another common sugar form
        Chem.MolFromSmarts("C1OC(CO)C(O)C1O"),    # Furanose form
        Chem.MolFromSmarts("C1OC(O)C(O)CC1"),     # Simple five-membered ring
    ]

    # Check for any sugar moiety
    if not any(mol.HasSubstructMatch(pattern) for pattern in sugar_patterns):
        return False, "No sugar moiety found"

    # Define generalized lipid pattern (a series of carbon atoms with some optional unsaturation)
    lipid_pattern = Chem.MolFromSmarts("[C,R0][C,R0][C,R0][C,R0][C,R0][C,R0][C]")  # Basic long aliphatic chain

    if not mol.HasSubstructMatch(lipid_pattern):
        return False, "No long carbon chain (lipid) found"

    return True, "Contains sugar moiety linked to a lipid component"