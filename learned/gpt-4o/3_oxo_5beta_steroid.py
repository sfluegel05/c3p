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

    # Update steroid backbone to a more general pattern that accounts for all rings
    steroid_pattern = Chem.MolFromSmarts("C12[C@@H]3C[C@@H]4[C@H]5CC[C@H](CCC5)[C@@]4(C)C[C@H]3CC2CCC1")  # Fused steroid rings
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found or wrong backbone configuration"
    
    # Check for specific 3-oxo group at correct position
    three_oxo_pattern = Chem.MolFromSmarts("C1=CC(=O)CC[C@H]1")  # More explicit 3-oxo pattern
    if not mol.HasSubstructMatch(three_oxo_pattern):
        return False, "No ketone group at 3rd position"

    # Confirm correct 5beta stereochemistry
    five_beta_pattern = Chem.MolFromSmarts("[C@@H](C)C(C)C")  # Updated to be more precise
    beta_matches = mol.GetSubstructMatches(five_beta_pattern)
    if not beta_matches or not any([match for match in beta_matches if mol.GetAtomWithIdx(match[0]).GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW]):
        return False, "5beta configuration not detected"

    return True, "Contains a 3-oxo group and 5beta configuration consistent with 3-oxo-5beta-steroid class"

# Example SMILES strings can be tested with this function for verification.