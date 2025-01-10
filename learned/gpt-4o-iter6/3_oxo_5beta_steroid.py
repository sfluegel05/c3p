"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
from rdkit import Chem

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5beta-steroid based on its SMILES string.
    A 3-oxo-5beta-steroid is any 3-oxo steroid that has beta-configuration at position 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-5beta-steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the 3-oxo group (carbonyl at position 3)
    oxo_pattern = Chem.MolFromSmarts("C(=O)[C@@H]")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No 3-oxo group found"

    # Check for 5beta configuration (beta hydrogen at position 5)
    beta_pattern = Chem.MolFromSmarts('[C@@H]1(CC[C@H](C1)[CH3])[CH3]')
    if not mol.HasSubstructMatch(beta_pattern):
        return False, "5beta configuration not found"

    # Look for steroid core structure (tetracyclic fused rings)
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4C3(CCC4)")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Tetracyclic steroid core not found"

    return True, "Contains 3-oxo group and 5beta-configuration in a steroid core"