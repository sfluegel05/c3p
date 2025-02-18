"""
Classifies: CHEBI:38131 lactol
"""
"""
Classifies: CHEBI_71648 lactol
"""
from rdkit import Chem

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol (cyclic hemiacetal) based on its SMILES string.
    A lactol is formed by intramolecular addition of a hydroxyl group to an aldehydic or ketonic carbonyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lactol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check all rings for the presence of a cyclic hemiacetal structure
    rings = mol.GetRingInfo().AtomRings()
    for ring in rings:
        ring_atoms = set(ring)
        # Check if the ring contains an oxygen
        if not any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring):
            continue
        
        # Iterate over each atom in the ring to find potential hemiacetal carbons
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:  # Focus on carbon atoms
                continue
            
            # Check if this carbon has a hydroxyl group
            has_oh = False
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() >= 1:
                    has_oh = True
                    break
            if not has_oh:
                continue
            
            # Check if this carbon is adjacent to a ring oxygen
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetIdx() in ring_atoms:
                    return True, "Cyclic hemiacetal structure detected"
    
    return False, "No lactol structure found"