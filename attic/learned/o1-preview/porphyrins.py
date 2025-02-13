"""
Classifies: CHEBI:26214 porphyrins
"""
"""
Classifies: CHEBI:8338 porphyrins
"""

from rdkit import Chem

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    A porphyrin contains a fundamental skeleton of four pyrrole nuclei united through the alpha-positions by four methine groups to form a macrocyclic structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a porphyrin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Flag to indicate if porphyrin ring is found
    porphyrin_found = False

    # Iterate over all rings
    for ring in atom_rings:
        # Check if ring is 16-membered
        if len(ring) == 16:
            # Count the number of nitrogen atoms in the ring
            num_nitrogens = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetSymbol() == 'N':
                    num_nitrogens += 1
            # Check if there are exactly 4 nitrogen atoms
            if num_nitrogens == 4:
                porphyrin_found = True
                break  # No need to check other rings

    if porphyrin_found:
        return True, "Porphyrin ring system detected (16-membered ring with 4 nitrogen atoms)"
    else:
        return False, "No porphyrin ring system detected"