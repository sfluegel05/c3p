"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
"""
Classifies: arenecarbaldehyde
Definition: Any aldehyde in which the carbonyl group is attached to an aromatic moiety.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    An arenecarbaldehyde has an aldehyde group (-CHO) directly attached to an aromatic ring.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_arenecarbaldehyde, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all aldehyde groups
    # [CX3H1](=O) matches a carbon with 3 connections, one H, and a double bond to oxygen
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    
    # Alternative pattern for explicit H representation
    aldehyde_pattern2 = Chem.MolFromSmarts("[CH1](=O)")
    
    matches = mol.GetSubstructMatches(aldehyde_pattern)
    matches2 = mol.GetSubstructMatches(aldehyde_pattern2)
    
    all_matches = list(set(matches + matches2))
    
    if not all_matches:
        return False, "No aldehyde group found"

    # For each aldehyde carbon, check if it's connected to an aromatic ring
    for match in all_matches:
        aldehyde_carbon = mol.GetAtomWithIdx(match[0])
        
        # Get all neighboring atoms of the aldehyde carbon
        for neighbor in aldehyde_carbon.GetNeighbors():
            # Skip the oxygen of the aldehyde group
            if neighbor.GetAtomicNum() == 8:
                continue
                
            # Check if the neighbor is aromatic
            if neighbor.GetIsAromatic():
                # Check if neighbor is part of an aromatic ring
                ring_info = mol.GetRingInfo()
                if ring_info.IsAtomInRing(neighbor.GetIdx()):
                    return True, "Contains aldehyde group directly attached to aromatic ring"
            
            # Check if neighbor is part of an aromatic ring system
            elif neighbor.GetAtomicNum() == 6:  # Carbon atom
                # Check if it's part of any ring
                ring_info = mol.GetRingInfo()
                if ring_info.IsAtomInRing(neighbor.GetIdx()):
                    # Check if any atoms in the same ring are aromatic
                    ring_atoms = set()
                    for ring in ring_info.AtomRings():
                        if neighbor.GetIdx() in ring:
                            ring_atoms.update(ring)
                    
                    for ring_atom_idx in ring_atoms:
                        if mol.GetAtomWithIdx(ring_atom_idx).GetIsAromatic():
                            return True, "Contains aldehyde group directly attached to aromatic ring system"

    return False, "No aldehyde group attached to aromatic ring"