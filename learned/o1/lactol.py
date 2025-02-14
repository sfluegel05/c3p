"""
Classifies: CHEBI:38131 lactol
"""
from rdkit import Chem

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol based on its SMILES string.
    A lactol is a cyclic hemiacetal formed by intramolecular addition of a hydroxy group to an aldehydic or ketonic carbonyl group, resulting in a ring containing an ether oxygen and a hydroxyl group on the adjacent carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lactol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ri = mol.GetRingInfo()
    atom_rings = ri.AtomRings()
    if not atom_rings:
        return False, "Molecule is not cyclic"

    found_lactol = False

    # Iterate over all rings
    for ring in atom_rings:
        ring_atoms = set(ring)
        # Find oxygen atoms in ring
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 8:
                # Oxygen in ring found
                # Get neighboring atoms
                neighbors = atom.GetNeighbors()
                for neighbor in neighbors:
                    if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() in ring_atoms:
                        # Neighbor is a carbon atom in the ring
                        # Check if this carbon has a hydroxyl group attached
                        has_OH = False
                        for nbr in neighbor.GetNeighbors():
                            if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in ring_atoms:
                                # Check if it's a hydroxyl oxygen (OH group)
                                if nbr.GetDegree() == 1:
                                    has_OH = True
                                    break
                        if has_OH:
                            return True, "Contains lactol functional group (cyclic hemiacetal)"
    return False, "Does not contain lactol functional group (cyclic hemiacetal)"