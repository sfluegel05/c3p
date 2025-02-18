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
    
    # SMARTS pattern for a cyclic hemiacetal: hydroxyl group attached to a carbon adjacent to ring oxygen
    # Matches [OH]-C-O in a ring structure
    lactol_pattern = Chem.MolFromSmarts('[OH]-C-O-*')
    if not mol.HasSubstructMatch(lactol_pattern):
        return False, "No hydroxyl adjacent to ring oxygen"
    
    # Verify that the oxygen is part of a ring
    rings = mol.GetRingInfo().AtomRings()
    for ring in rings:
        ring_atoms = set(ring)
        # Check for oxygen in the ring
        oxygen_in_ring = any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring)
        if not oxygen_in_ring:
            continue
        # Check for a carbon in the ring with hydroxyl group adjacent to the ring oxygen
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:  # Oxygen atom in the ring
                # Get neighboring atoms in the ring
                neighbors = [n for n in atom.GetNeighbors() if n.GetIdx() in ring_atoms]
                for neighbor in neighbors:
                    if neighbor.GetAtomicNum() == 6:  # Carbon adjacent to ring oxygen
                        # Check if this carbon has a hydroxyl group
                        for bond in neighbor.GetBonds():
                            other = bond.GetOtherAtom(neighbor)
                            if other.GetAtomicNum() == 8 and other.GetTotalNumHs() >= 1:
                                return True, "Cyclic hemiacetal structure detected"
    
    return False, "No lactol structure found"