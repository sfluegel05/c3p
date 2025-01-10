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

    # Define a basic pattern for sugar moiety - usually pyranose (6-member) or furanose (5-member) rings
    sugar_pattern = Chem.MolFromSmarts("C1(O[CH](CO)[C@@H](O)[C@H](O)[C@H]1O) | C1(O[CH](CO)[C@@H](O)[C@H](O)[C@H](O)[C@@H]1O)")

    # Check for sugar moiety
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar moiety found"

    # Check for long carbon chains representing lipids
    lipid_pattern = Chem.MolFromSmarts("CCCCCCCCCCC")  # Simple pattern for a long carbon chain

    if not mol.HasSubstructMatch(lipid_pattern):
        return False, "No long carbon chain (lipid) found"

    return True, "Contains sugar moiety linked to a lipid component"