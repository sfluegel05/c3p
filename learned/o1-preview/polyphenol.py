"""
Classifies: CHEBI:26195 polyphenol
"""
"""
Classifies: polyphenol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polyphenol(smiles: str):
    """
    Determines if a molecule is a polyphenol based on its SMILES string.
    A polyphenol is a molecule containing 2 or more benzene rings each substituted by at least one hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyphenol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure aromaticity perception
    try:
        Chem.SanitizeMol(mol)
    except Chem.SanitizeException as e:
        return False, f"Sanitization failed: {e}"

    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    benzene_rings = []
    for ring in atom_rings:
        if len(ring) != 6:
            continue
        is_benzene = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                is_benzene = False
                break
            if not atom.GetIsAromatic():
                is_benzene = False
                break
        if is_benzene:
            benzene_rings.append(ring)

    if len(benzene_rings) == 0:
        return False, "No benzene rings found"

    # Keep track of benzene rings substituted by at least one hydroxy group
    substituted_rings = 0

    for ring in benzene_rings:
        ring_atom_indices = set(ring)
        has_hydroxy = False

        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                # Skip atoms within the ring
                if neighbor_idx in ring_atom_indices:
                    continue
                # Check if neighbor is an oxygen connected via single bond (hydroxy group)
                if neighbor.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(atom_idx, neighbor_idx).GetBondType() == Chem.rdchem.BondType.SINGLE:
                    # Check if oxygen has a hydrogen (i.e., is a hydroxy group)
                    if neighbor.GetTotalNumHs(includeNeighbors=True) > 0:
                        has_hydroxy = True
                        break
            if has_hydroxy:
                break

        if has_hydroxy:
            substituted_rings += 1

    if substituted_rings >= 2:
        return True, f"Contains {substituted_rings} benzene rings each substituted with at least one hydroxy group"
    else:
        return False, f"Only {substituted_rings} benzene rings substituted with hydroxy groups"

__metadata__ = {
    'chemical_class': {
        'name': 'polyphenol',
        'definition': 'Members of the class of phenols that contain 2 or more benzene rings each of which is substituted by at least one hydroxy group.'
    }
}