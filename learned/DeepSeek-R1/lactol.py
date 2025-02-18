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
    
    # Check all rings for lactol pattern: adjacent ring oxygen and hydroxyl-bearing carbon
    rings = mol.GetRingInfo().AtomRings()
    for ring in rings:
        ring_atoms = set(ring)
        for i in range(len(ring)):
            atom_idx1 = ring[i]
            atom_idx2 = ring[(i+1) % len(ring)]  # Next atom in the ring
            
            atom1 = mol.GetAtomWithIdx(atom_idx1)
            atom2 = mol.GetAtomWithIdx(atom_idx2)
            
            # Check for O-C(OH) adjacency in the ring
            if (atom1.GetAtomicNum() == 8 and atom2.GetAtomicNum() == 6):
                # Check if atom2 has exactly one hydroxyl group and total two oxygen neighbors
                hydroxyls = 0
                oxygen_neighbors = 0
                for nbr in atom2.GetNeighbors():
                    if nbr.GetAtomicNum() == 8:
                        oxygen_neighbors += 1
                        if nbr.GetTotalNumHs() >= 1:
                            hydroxyls += 1
                if hydroxyls == 1 and oxygen_neighbors == 2:
                    return True, "Adjacent ring oxygen and hydroxyl-bearing carbon (lactol)"
            elif (atom2.GetAtomicNum() == 8 and atom1.GetAtomicNum() == 6):
                # Check atom1 similarly
                hydroxyls = 0
                oxygen_neighbors = 0
                for nbr in atom1.GetNeighbors():
                    if nbr.GetAtomicNum() == 8:
                        oxygen_neighbors += 1
                        if nbr.GetTotalNumHs() >= 1:
                            hydroxyls += 1
                if hydroxyls == 1 and oxygen_neighbors == 2:
                    return True, "Adjacent ring oxygen and hydroxyl-bearing carbon (lactol)"
    
    return False, "No lactol structure detected"