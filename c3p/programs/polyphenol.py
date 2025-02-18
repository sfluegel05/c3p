"""
Classifies: CHEBI:26195 polyphenol
"""
"""
Classifies: Polyphenol
"""
from rdkit import Chem

def is_polyphenol(smiles: str):
    """
    Determines if a molecule is a polyphenol based on its SMILES string.
    A polyphenol is a member of the class of phenols that contains two or more benzene rings,
    each of which is substituted by at least one hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a polyphenol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure proper aromaticity perception
    Chem.SanitizeMol(mol)

    # Define phenol SMARTS pattern (hydroxy group attached to aromatic carbon)
    phenol_pattern = Chem.MolFromSmarts('[OX2H]-c')

    # Find all phenol groups in the molecule
    phenol_matches = mol.GetSubstructMatches(phenol_pattern)

    if not phenol_matches:
        return False, "No phenol groups found"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Set to keep track of unique benzene rings with hydroxy substitution
    benzene_rings_with_hydroxy = set()

    # Iterate over phenol groups
    for match in phenol_matches:
        oxygen_idx, carbon_idx = match  # Indices of O and C in [OH]-c

        # Get the rings that the carbon atom is part of
        rings_with_carbon = [ring for ring in atom_rings if carbon_idx in ring]

        # Check each ring
        for ring in rings_with_carbon:
            # Check if ring is six-membered
            if len(ring) != 6:
                continue

            # Check if all atoms in the ring are aromatic
            is_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
            if not is_aromatic:
                continue  # Not an aromatic ring

            # Add the ring (as a frozenset of atom indices) to the set
            benzene_rings_with_hydroxy.add(frozenset(ring))

    # Determine the number of unique benzene rings with hydroxy substitution
    ring_count = len(benzene_rings_with_hydroxy)

    if ring_count >= 2:
        return True, f"Molecule is a polyphenol with {ring_count} benzene rings each substituted by at least one hydroxy group"
    else:
        return False, f"Molecule has {ring_count} benzene rings with hydroxy substitution, needs at least 2 to be a polyphenol"