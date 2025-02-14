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
    
    This function checks if the molecule contains an aromatic ring (benzene ring)
    with at least four chlorine substituents attached to the ring atoms.

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

    found = False  # Flag indicating if tetrachlorobenzene ring is found

    # Get ring information
    ri = mol.GetRingInfo()
    # Get list of all rings (list of atom indices)
    ring_list = ri.AtomRings()

    # Iterate over all rings
    for ring in ring_list:
        # Check if ring is aromatic
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if not all(atom.GetIsAromatic() for atom in ring_atoms):
            continue
        chlorine_count = 0  # Number of chlorines attached to the ring
        # Examine substituents on ring atoms
        for atom in ring_atoms:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in ring:
                    continue  # Skip ring atoms
                if neighbor.GetAtomicNum() == 17:
                    chlorine_count += 1  # Chlorine substituent
        # Check if there are at least 4 chlorines attached to the ring
        if chlorine_count >= 4:
            found = True
            break  # Found tetrachlorobenzene ring

    if found:
        return True, "Contains aromatic ring with at least four chlorine substituents"
    else:
        return False, "Does not contain aromatic ring with at least four chlorine substituents"