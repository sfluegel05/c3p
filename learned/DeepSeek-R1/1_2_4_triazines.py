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
    at positions 1, 2, and 4 of the core ring structure.
    
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
    
    # Get all rings in the molecule
    rings = mol.GetRingInfo().AtomRings()
    
    for ring in rings:
        if len(ring) != 6:
            continue  # Only consider six-membered rings
        
        # Check if all atoms in the ring are aromatic
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
            
        # Get nitrogen positions in the ring (0-5 index)
        n_positions = [i for i, idx in enumerate(ring) 
                      if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7]
        
        if len(n_positions) != 3:
            continue  # Need exactly three nitrogens
        
        # Check for 1,2,4 pattern: two adjacent nitrogens and third two positions away
        for i in range(len(n_positions)):
            for j in range(i+1, len(n_positions)):
                a, b = sorted([n_positions[i], n_positions[j]])
                
                # Check adjacency (direct or wrap-around)
                if (b - a == 1) or (a == 0 and b == 5):
                    # Calculate expected third nitrogen position
                    third = (b + 2) % 6
                    if third in n_positions:
                        return True, "Contains six-membered aromatic ring with three nitrogens in 1,2,4 positions"
    
    return False, "No 1,2,4-triazine skeleton found"