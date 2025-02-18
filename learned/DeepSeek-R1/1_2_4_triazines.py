"""
Classifies: CHEBI:39410 1,2,4-triazines
"""
"""
Classifies: 1,2,4-triazines
"""
from rdkit import Chem

def is_1_2_4_triazines(smiles: str):
    """
    Determines if a molecule contains a 1,2,4-triazine skeleton.
    A 1,2,4-triazine has a six-membered aromatic ring with three nitrogen atoms
    arranged such that two are adjacent and the third is two positions away.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains 1,2,4-triazine skeleton, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    try:
        Chem.SanitizeMol(mol)
    except:
        return False, "Sanitization failed"
    
    # Get all rings and check each six-membered aromatic ring
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) != 6:
            continue
        
        # Check aromaticity for all atoms in the ring
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        
        # Collect nitrogen positions in the ring's atom order
        n_positions = [i for i, atom_idx in enumerate(ring) 
                      if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 7]
        
        if len(n_positions) != 3:
            continue  # Need exactly three nitrogens
        
        # Check all combinations of nitrogen pairs
        for i in range(len(n_positions)):
            for j in range(i+1, len(n_positions)):
                a, b = sorted([n_positions[i], n_positions[j]])
                
                # Check adjacency (direct or wrap-around)
                if (b - a == 1) or (a == 0 and b == 5):
                    # Check third nitrogen two positions after b
                    third = (b + 2) % 6
                    if third in n_positions:
                        return True, "Contains 1,2,4-triazine skeleton"
                    
                    # Check third nitrogen two positions before a (reverse direction)
                    third_rev = (a - 2) % 6
                    if third_rev in n_positions:
                        return True, "Contains 1,2,4-triazine skeleton"
    
    return False, "No 1,2,4-triazine skeleton found"