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

    # Update to a more general steroid backbone pattern
    steroid_pattern = Chem.MolFromSmarts("[C@]12CC[C@@H]3[C@@H](CC[C@@H]4CC(=O)CC[C@]43C)[C@@H]1CC[C@@H]2O")  # Generic 5-beta steroid structure
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found or wrong backbone configuration"
    
    # Check for 3-oxo group at A-ring
    three_oxo_pattern = Chem.MolFromSmarts("C(=O)[C@]1([H])CC")  # Captures a typical A-ring with 3-oxo
    if not mol.HasSubstructMatch(three_oxo_pattern):
        return False, "No ketone group at 3rd position within steroid core"

    # Confirm correct 5beta stereochemistry in the B-ring context
    five_beta_pattern = Chem.MolFromSmarts("C12CCC[C@@H]1")
    five_beta_matches = mol.GetSubstructMatches(five_beta_pattern)
    if not five_beta_matches:
        return False, "5beta configuration not detected"

    return True, "Contains a 3-oxo group and 5beta configuration consistent with 3-oxo-5beta-steroid class"