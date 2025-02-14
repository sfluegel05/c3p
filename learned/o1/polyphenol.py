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

    # Initialize count of benzene rings with hydroxy substitutions
    benzene_with_hydroxy_count = 0

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Iterate over all rings in the molecule
    for ring in atom_rings:
        # Check if ring is six-membered
        if len(ring) != 6:
            continue

        # Check if all atoms in the ring are aromatic (benzene ring)
        is_benzene = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if not atom.GetIsAromatic():
                is_benzene = False
                break
        if not is_benzene:
            continue  # Not a benzene ring

        # Check for hydroxy substitution on the benzene ring
        has_hydroxy = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for neighbor in atom.GetNeighbors():
                nbr_idx = neighbor.GetIdx()
                # Skip if neighbor is part of the ring
                if nbr_idx in ring:
                    continue
                # Check if neighbor is oxygen
                if neighbor.GetAtomicNum() == 8:
                    # Check if the oxygen has only one heavy atom neighbor (the ring atom)
                    if neighbor.GetDegree() == 1:
                        has_hydroxy = True
                        break  # Found hydroxy group
            if has_hydroxy:
                break  # No need to check other atoms in the ring

        if has_hydroxy:
            benzene_with_hydroxy_count += 1

    # Determine if molecule is a polyphenol
    if benzene_with_hydroxy_count >= 2:
        return True, f"Molecule is a polyphenol with {benzene_with_hydroxy_count} benzene rings each substituted by at least one hydroxy group"
    else:
        return False, f"Molecule has {benzene_with_hydroxy_count} benzene rings with hydroxy substitution, needs at least 2 to be a polyphenol"