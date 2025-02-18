"""
Classifies: CHEBI:26195 polyphenol
"""
"""
Classifies: CHEBI:26195 polyphenol
"""
from rdkit import Chem

def is_polyphenol(smiles: str):
    """
    Determines if a molecule is a polyphenol based on its SMILES string.
    A polyphenol has at least two benzene rings, each substituted with at least one hydroxyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyphenol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all benzene rings (6-membered aromatic rings)
    benzene_rings = []
    rings = mol.GetRingInfo().AtomRings()
    for ring in rings:
        if len(ring) == 6:
            is_aromatic = all(mol.GetAtomWithIdx(atom).GetIsAromatic() for atom in ring)
            if is_aromatic:
                benzene_rings.append(ring)

    if len(benzene_rings) < 2:
        return False, f"Found {len(benzene_rings)} benzene rings, need at least 2"

    # Find hydroxyl groups attached to aromatic carbons
    hydroxyl_pattern = Chem.MolFromSmarts('[c]-[OH1]')  # Matches -OH on aromatic carbons
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    hydroxyl_carbons = {match[0] for match in hydroxyl_matches}

    # Check each benzene ring for at least one hydroxyl substituent
    valid_rings = 0
    for ring in benzene_rings:
        for atom in ring:
            if atom in hydroxyl_carbons:
                valid_rings += 1
                break  # Only one hydroxyl needed per ring

    if valid_rings >= 2:
        return True, f"Found {valid_rings} benzene rings with hydroxyl groups"
    else:
        return False, f"Only {valid_rings} benzene rings with hydroxyl groups"