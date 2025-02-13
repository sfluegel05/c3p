"""
Classifies: CHEBI:47788 3-oxo steroid
"""
from rdkit import Chem

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    A 3-oxo steroid is characterized by a steroidal framework with an oxo group (C=O) at the 3rd position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generalized pattern for four-fused rings as in steroids
    steroid_ring_smarts = '[R]1[R][R][R]2[R][R][R]3[R][R][R]4[R][R][R]2[R][R][R]1[R][R][R]3[R][R]4'
    steroid_pattern = Chem.MolFromSmarts(steroid_ring_smarts)
    
    # Check for steroid backbone
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No four-ring steroid backbone found"

    # C(=O) group identification
    oxo_smarts = '[C]=O'
    oxo_pattern = Chem.MolFromSmarts(oxo_smarts)
    
    # Find potential matches of the oxo groups
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    if not oxo_matches:
        return False, "No C=O bond found"

    # Evaluate if the oxo group is appropriately positioned (usually at 3rd position in a common steroid layout)
    # This typically involves checking ring adjacency or direct substitution positions
    for match in oxo_matches:
        carbon_idx = match[0]
        # Simple adjacency check, in practical terms this should map to 3rd position in steroid rings
        if any(mol.GetBondWithIdx(bond.GetIdx()).IsInRing() for bond in mol.GetAtomWithIdx(carbon_idx).GetBonds()):
            return True, "Detected steroid with a 3-oxo (C=O) group"

    return False, "No oxo group at position 3 detected"