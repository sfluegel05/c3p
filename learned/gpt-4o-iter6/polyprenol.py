"""
Classifies: CHEBI:26199 polyprenol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    A polyprenol is a prenol composed of more than one isoprene unit, ending with a hydroxyl group.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for isoprene units including stereochemistry
    isoprene_pattern = Chem.MolFromSmarts("[C:1](=[C:2])[C:3](C)[CH2]")
    
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 2:
        return False, f"Found {len(isoprene_matches)} isoprene units, need more than 1"
    
    # Identify hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    
    if not hydroxyl_matches:
        return False, "No hydroxyl group found"
    
    # Check connectivity of hydroxyl group
    for match in hydroxyl_matches:
        oh_atom_idx = match[0]
        oh_atom = mol.GetAtomWithIdx(oh_atom_idx)
        neighbor_indices = [neighbor.GetIdx() for neighbor in oh_atom.GetNeighbors()]

        # Check if the hydroxyl is at the end of isoprene unit chains
        for index in neighbor_indices:
            if any(index in m for m in isoprene_matches):
                # Detected OH connected reliably at the chain end
                return True, "Contains more than one isoprene unit with hydroxyl end"
    
    return False, "The hydroxyl group is not positioned correctly relative to isoprene units"