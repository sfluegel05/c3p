"""
Classifies: CHEBI:28966 chlorophyll
"""
"""
Classifies: chlorophyll

Definition: A family of magnesium porphyrins, defined by the presence of a fifth ring beyond the four pyrrole-like rings. The rings can have various side chains which usually include a long phytol chain.
"""
from rdkit import Chem

def is_chlorophyll(smiles: str):
    """
    Determines if a molecule is a chlorophyll based on its SMILES string.
    A chlorophyll is characterized by a magnesium porphyrin core with a fifth ring and side chains,
    often including a long phytol chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chlorophyll, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for magnesium atom (atomic number 12)
    mg_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 12]
    if not mg_atoms:
        return False, "No magnesium atom found"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Find macrocyclic rings (size >=16) containing at least 4 nitrogen atoms
    macrocycle_found = False
    for ring in atom_rings:
        if len(ring) >= 16:
            # Check if ring contains at least 4 nitrogen atoms
            nitrogen_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            if nitrogen_count >= 4:
                macrocycle_found = True
                macrocycle_ring = ring
                break
    if not macrocycle_found:
        return False, "No macrocyclic ring with at least 4 nitrogen atoms found"

    # Check that magnesium is connected to nitrogen atoms in the macrocycle
    macrocycle_atom_indices = set(macrocycle_ring)
    for mg in mg_atoms:
        neighbors = mg.GetNeighbors()
        # Find nitrogen neighbors that are part of the macrocycle
        nitrogen_neighbors = [atom for atom in neighbors if atom.GetAtomicNum() == 7 and atom.GetIdx() in macrocycle_atom_indices]
        if len(nitrogen_neighbors) >= 4:
            # Passed all checks
            return True, "Contains magnesium porphyrin core with fifth ring"
    return False, "Magnesium not coordinated to nitrogen atoms in macrocycle"

__metadata__ = {
    'chemical_class': {
        'name': 'chlorophyll',
        'definition': 'A family of magnesium porphyrins, defined by the presence of a fifth ring beyond the four pyrrole-like rings. The rings can have various side chains which usually include a long phytol chain.'
    }
}