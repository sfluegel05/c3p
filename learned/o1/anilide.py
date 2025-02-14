"""
Classifies: CHEBI:13248 anilide
"""
"""
Classifies: anilide
"""
from rdkit import Chem

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.
    An anilide is any aromatic amide obtained by acylation of aniline,
    meaning it has an amide group where the nitrogen is directly attached
    to a benzene ring (aromatic ring of six carbons).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anilide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify amide functional groups: [NX3][CX3](=O)
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    
    if not amide_matches:
        return False, "No amide functional group found"
    
    ring_info = mol.GetRingInfo()
    for match in amide_matches:
        nitrogen_idx, carbonyl_c_idx = match[0], match[1]
        nitrogen = mol.GetAtomWithIdx(nitrogen_idx)
        carbonyl_c = mol.GetAtomWithIdx(carbonyl_c_idx)
        
        # Check if nitrogen is not in a ring
        if ring_info.IsAtomInRingOfSize(nitrogen_idx, 5) or ring_info.IsAtomInRingOfSize(nitrogen_idx, 6):
            continue  # Skip if nitrogen is in a ring
        
        # Find atoms connected to the nitrogen, excluding the carbonyl carbon
        neighbor_atoms = [a for a in nitrogen.GetNeighbors() if a.GetIdx() != carbonyl_c_idx]
        
        for neighbor in neighbor_atoms:
            # Check if neighbor is an aromatic carbon
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIsAromatic():
                neighbor_idx = neighbor.GetIdx()
                # Check if the aromatic carbon is part of a benzene ring
                is_benzene = False
                for ring in ring_info.AtomRings():
                    if neighbor_idx in ring and len(ring) == 6:
                        # Check if all atoms in the ring are carbons
                        if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
                            is_benzene = True
                            break
                if is_benzene:
                    return True, "Contains an anilide moiety (amide nitrogen attached to benzene ring)"
    return False, "Does not contain an anilide moiety"