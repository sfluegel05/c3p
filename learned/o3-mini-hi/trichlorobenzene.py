"""
Classifies: CHEBI:27096 trichlorobenzene
"""
"""
Classifies: trichlorobenzene
Defined as any member of the class of chlorobenzenes carrying three chloro substituents at unspecified positions.
"""

from rdkit import Chem

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule is a trichlorobenzene based on its SMILES string.
    A trichlorobenzene is defined as a molecule that contains a benzene ring bearing exactly
    three chlorine atoms as peripheral substituents.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a trichlorobenzene, False otherwise
        str: A reason explaining the classification decision.
    """
    # Parse the SMILES string to a molecule object using RDKit.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Retrieve information about rings in the molecule.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Loop over each ring found.
    for ring in atom_rings:
        # We are interested only in six-membered rings.
        if len(ring) != 6:
            continue

        # Check if all atoms in the ring are aromatic (benzene-like)
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue

        # For a valid aromatic benzene ring, count the number of chlorine substituents.
        chloro_count = 0
        # Iterate through each atom in the ring.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Examine all neighbors of the ring atom which are not part of the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() not in ring and nbr.GetAtomicNum() == 17:  # Chlorine atomic number is 17
                    chloro_count += 1

        # If this aromatic ring has exactly three Cl substituents, we classify as trichlorobenzene.
        if chloro_count == 3:
            return True, "Found a benzene ring with exactly three chloro substituents"

    return False, "No benzene ring with exactly three chloro substituents found"