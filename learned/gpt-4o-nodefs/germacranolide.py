"""
Classifies: CHEBI:73011 germacranolide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_germacranolide(smiles: str):
    """
    Determines if a molecule is a germacranolide based on its SMILES string.
    Germacranolides are sesquiterpene lactones with a cyclodecene framework
    and a lactone group.

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
    
    # Check for a 10-membered ring
    ring_info = mol.GetRingInfo()
    found_10_membered_ring = False
    for ring in ring_info.AtomRings():
        if len(ring) == 10:
            found_10_membered_ring = True
            break

    if not found_10_membered_ring:
        return False, "No 10-membered ring found, typical for germacranolide"

    # Look for a lactone group (cyclic ester with at least one =O in C=O)
    lactone_pattern = Chem.MolFromSmarts("O=C1OC=CC=CC=CC(C)=CC1")
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "Lactone group not found"

    # Check for multiple double bonds
    num_double_bonds = rdMolDescriptors.CalcNumDoubleBonds(mol)
    if num_double_bonds < 2:
        return False, "Insufficient double bonds, typical germacranolides are polyunsaturated"

    # Assess presence of oxygen functionalities: hydroxyl and ester groups
    oxy_patterns = [Chem.MolFromSmarts("C=O"), Chem.MolFromSmarts("O[C@H]"), Chem.MolFromSmarts("O[C@@H]")]
    for pattern in oxy_patterns:
        if not mol.HasSubstructMatch(pattern):
            return False, "Unmatched oxygen functionalities"

    # Passed all checks: likely a germacranolide
    return True, "Matches typical structure of a germacranolide"