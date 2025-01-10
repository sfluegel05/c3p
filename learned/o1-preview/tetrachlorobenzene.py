"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
"""
Classifies: tetrachlorobenzene
"""
from rdkit import Chem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule is a tetrachlorobenzene based on its SMILES string.
    A tetrachlorobenzene is any chlorobenzene carrying four chloro groups at unspecified positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrachlorobenzene, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ri = mol.GetRingInfo()
    found = False  # Flag indicating if tetrachlorobenzene ring is found
    
    # Iterate over all atom rings
    for ring in ri.AtomRings():
        # Check if ring is size 6 and aromatic
        if len(ring) != 6:
            continue
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if not all(atom.GetAtomicNum() == 6 and atom.GetIsAromatic() for atom in ring_atoms):
            continue
        # Check that ring is not fused (each atom is only in one ring)
        if any(ri.NumAtomRings(idx) > 1 for idx in ring):
            continue
        chlorine_count = 0
        other_substituents = False
        # Examine substituents on ring atoms
        for atom in ring_atoms:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in ring:
                    continue  # Skip ring atoms
                atomic_num = neighbor.GetAtomicNum()
                if atomic_num == 17:
                    chlorine_count += 1  # Chlorine substituent
                elif atomic_num == 1:
                    continue  # Hydrogen, ignore
                else:
                    # Found substituent other than chlorine or hydrogen
                    other_substituents = True
                    break
            if other_substituents:
                break
        if other_substituents:
            continue  # Skip to next ring
        # Check if there are exactly 4 chlorines attached to the ring
        if chlorine_count == 4:
            found = True
            break  # Found tetrachlorobenzene ring
    
    if found:
        return True, "Contains isolated benzene ring with exactly four chlorine substituents"
    else:
        return False, "Does not contain isolated benzene ring with exactly four chlorine substituents"