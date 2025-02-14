"""
Classifies: CHEBI:26660 sesterterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    A sesterterpenoid is derived from sesterterpenes composed of five isoprene units (C25 skeleton),
    possibly modified by rearrangement or removal of small groups (e.g., methyl groups).
    They typically contain ring structures, methyl branching, and oxygen-containing functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesterterpenoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20 or c_count > 30:
        return False, f"Carbon count is {c_count}, not typical for sesterterpenoids (expected ~25 carbons)"
    
    # Check for ring structures
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings == 0:
        return False, "No ring structures found; sesterterpenoids often contain rings"
    
    # Check for methyl branching (tertiary and quaternary carbons)
    tertiary_carbons = 0
    quaternary_carbons = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            neighbor_carbons = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 6)
            if neighbor_carbons == 3:
                tertiary_carbons += 1
            elif neighbor_carbons == 4:
                quaternary_carbons += 1
    if (tertiary_carbons + quaternary_carbons) < 3:
        return False, f"Found {tertiary_carbons} tertiary and {quaternary_carbons} quaternary carbons; sesterterpenoids typically have methyl branching"
    
    # Check for oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count == 0:
        return False, "No oxygen atoms found; sesterterpenoids often contain oxygenated functional groups"
    
    # Check for typical terpenoid functional groups
    # e.g., hydroxyl, carbonyl, ester, ether
    functional_groups = [
        Chem.MolFromSmarts('[OX2H]'),     # Hydroxyl group
        Chem.MolFromSmarts('C=O'),        # Carbonyl group
        Chem.MolFromSmarts('[CX3](=O)[OX2H1]'),  # Carboxylic acid
        Chem.MolFromSmarts('[CX3](=O)[OX2][CX3]'), # Ester
        Chem.MolFromSmarts('[OX2][CX4][OX2]'),    # Ether
    ]
    fg_found = False
    for fg in functional_groups:
        if mol.HasSubstructMatch(fg):
            fg_found = True
            break
    if not fg_found:
        return False, "No typical terpenoid functional groups found"

    # Check for isoprene units (C5H8)
    # It's challenging to identify isoprene units directly, but we can look for substructures
    # We can check for units with patterns like C=C-C-C=C
    isoprene_pattern = Chem.MolFromSmarts('C=C-C-C=C')
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 1:
        return False, "No isoprene-like units found"

    return True, "Molecule meets criteria for a sesterterpenoid"