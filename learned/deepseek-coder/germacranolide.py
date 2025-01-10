"""
Classifies: CHEBI:73011 germacranolide
"""
"""
Classifies: CHEBI:37698 germacranolide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_germacranolide(smiles: str):
    """
    Determines if a molecule is a germacranolide based on its SMILES string.
    A germacranolide is a sesquiterpene lactone based on the germacrane skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a germacranolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for 10-membered ring (germacrane skeleton)
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    if 10 not in ring_sizes:
        return False, "No 10-membered ring (germacrane skeleton) found"

    # Look for lactone group (cyclic ester) in a more specific configuration
    lactone_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CX4]1[CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4]1")
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone group in typical germacranolide configuration found"

    # Check for sesquiterpene base (approximately 15 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15 or c_count > 25:  # Allow more flexibility for functional groups
        return False, f"Not a sesquiterpene (15-25 carbons expected, found {c_count})"

    # Check for typical germacranolide connectivity
    # Look for the characteristic fused ring system with lactone
    germacranolide_pattern = Chem.MolFromSmarts("[CX4]1[CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4]1.[CX3](=[OX1])[OX2][CX4]")
    if not mol.HasSubstructMatch(germacranolide_pattern):
        return False, "Does not match typical germacranolide connectivity pattern"

    return True, "Contains 10-membered ring, lactone group, and sesquiterpene base structure"