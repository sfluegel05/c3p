"""
Classifies: CHEBI:26935 tetraterpenoid
"""
from rdkit import Chem

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    A tetraterpenoid typically contains extensive conjugated double bonds,
    possibly cyclic structures, often exceeding 40 carbon atoms, and diverse
    functional groups.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is a tetraterpenoid, False otherwise
        str: Reason for classification or rejection
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check approximate carbon count
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 30 or c_count > 60:
        return False, f"Unsupported carbon count: {c_count} (expected 30-60)"
    
    # Look for conjugated double bond systems
    conjugated_bond_pattern = Chem.MolFromSmarts("C=C-C=C-C=C")  # Flexible pattern
    if not mol.HasSubstructMatch(conjugated_bond_pattern):
        return False, "Insufficient conjugated double-bond systems"
    
    # Check for typical functional groups like ketones, hydroxyls, ethers
    functional_group_detected = False
    functional_group_smarts = [
        "[CX3]=[OX1]",  # C=O for ketones
        "[OH]",         # Hydroxyl groups
        "[CX4][OX2][CX4]"  # Ethers
    ]
    
    for fg_smarts in functional_group_smarts:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(fg_smarts)):
            functional_group_detected = True
            break
    
    if not functional_group_detected:
        return False, "Lacking typical functional groups found in tetraterpenoids"
    
    # Check for possible cyclic structures
    ring_count = mol.GetRingInfo().NumRings()
    # Allow flexible presence or absence of rings
    ring_criterion_met = (ring_count > 0) or (ring_count == 0)
    
    if not ring_criterion_met:
        return False, "No cyclic structures detected, important for tetraterpenoids"
    
    return True, "Meets criteria for a tetraterpenoid based on structural attributes and functional group presence"

# Example usage
example_smiles = "O=C1C(=C(/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/[C@@H](OC)C(O)(C)C)\\C)\\C)/C)/C)C(C)(C)CC1)C"
print(is_tetraterpenoid(example_smiles))