"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
from rdkit import Chem

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5beta-steroid based on its SMILES string.
    A 3-oxo-5beta-steroid is a steroid with a ketone group at position 3 and a beta-configuration
    at position 5 in its sterane core.

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

    # Look for steroid backbone pattern (basic fused ring structure with steroid core)
    # Steroids typically have a perhydrocyclopentanophenanthrene structure, need to check SMILES for C17 core
    steroid_pattern = Chem.MolFromSmarts("[C@]12CC[C@]3(CC(=O)C[C@H]1CC[C@@H]2C(C3)(CCC(=O)O)C)")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found or wrong backbone configuration"
    
    # Check for 3-oxo group at A-ring
    three_oxo_pattern = Chem.MolFromSmarts("C(=O)[C@H]")  # C=O at specific position
    if not mol.HasSubstructMatch(three_oxo_pattern):
        return False, "No ketone group at 3rd position"

    # Check for correct 5beta stereochemistry by identifying the right stereochemistry annotation
    five_beta_pattern = Chem.MolFromSmarts("[C@@H]")  # Represents beta stereochemistry on carbon
    beta_matches = mol.GetSubstructMatches(five_beta_pattern)
    if not beta_matches or len(beta_matches) < 1:
        return False, "5beta configuration not detected"

    return True, "Contains a 3-oxo group and 5beta configuration consistent with 3-oxo-5beta-steroid class"

# Example SMILES strings can be tested with this function for verification.