"""
Classifies: CHEBI:26935 tetraterpenoid
"""
from rdkit import Chem

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    A tetraterpenoid typically has conjugated double bonds, cyclic structures, 
    and certain functional groups, usually spanning long carbon chains.
    
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
    
    # Check approximate carbon count, but allow a wider range
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 35 or c_count > 55:
        return False, f"Unsupported carbon count: {c_count} (expected 35-55)"
    
    # Look for extended conjugated double bond systems
    polyene_pattern = Chem.MolFromSmarts("C=C" * 10)  # Look for long polyene chains
    if not mol.HasSubstructMatch(polyene_pattern):
        return False, "Insufficient or not extensive conjugated double-bond systems"
    
    # Check for typical functional groups like ketones and hydroxyls
    functional_group_detected = False
    functional_group_types = ["C=O", "[OH]"]
    
    for fg in functional_group_types:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(fg)):
            functional_group_detected = True
    
    if not functional_group_detected:
        return False, "Lacking typical functional groups found in tetraterpenoids"
    
    # Check for cyclic structures often present in tetraterpenoids
    ring_count = mol.GetRingInfo().NumRings()
    if ring_count == 0:
        return False, "No cyclic structures detected, important for tetraterpenoids"
    
    return True, "Meets criteria for a tetraterpenoid based on structural attributes and functional group presence"

# Example usage
example_smiles = "O=C1C(=C(/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/[C@@H](OC)C(O)(C)C)\\C)\\C)/C)/C)C(C)(C)CC1)C"
print(is_tetraterpenoid(example_smiles))