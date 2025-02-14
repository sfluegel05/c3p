"""
Classifies: CHEBI:47916 flavonoid
"""
"""
Classifies: Flavonoid
"""

from rdkit import Chem

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    A flavonoid is based on 1-benzopyran (chromene) with an aryl substituent at position 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define chromene core SMARTS pattern (1-benzopyran)
    chromene_smarts = 'c1ccc2occ(c2c1)'  # Chromene core
    chromene_pattern = Chem.MolFromSmarts(chromene_smarts)
    if chromene_pattern is None:
        return False, "Failed to construct chromene SMARTS pattern"

    # Find matches of chromene core in molecule
    chromene_matches = mol.GetSubstructMatches(chromene_pattern)
    if not chromene_matches:
        return False, "Molecule does not contain chromene core"

    # For each match of chromene core, check for aryl substituent at position 2
    # In the chromene SMARTS, the atom at index 3 corresponds to position 2
    for match in chromene_matches:
        pos2_atom_idx = match[3]  # Atom index at position 2 of chromene
        pos2_atom = mol.GetAtomWithIdx(pos2_atom_idx)

        # Check if the atom at position 2 has a neighbor not in the chromene core
        for neighbor in pos2_atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx in match:
                continue  # Skip atoms that are part of the chromene core

            # Check if the neighbor is part of an aromatic ring (aryl group)
            is_aryl = False
            for ring in neighbor.GetOwningMol().GetRingInfo().AtomRings():
                if neighbor_idx in ring:
                    # Check if the ring is aromatic and six-membered
                    if len(ring) == 6 and all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                        is_aryl = True
                        break
            if is_aryl:
                return True, "Molecule contains flavonoid core with aryl substituent at position 2"

    return False, "Molecule does not contain aryl substituent at position 2 of chromene core"