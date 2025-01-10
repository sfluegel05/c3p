"""
Classifies: CHEBI:35915 sterol ester
"""
"""
Classifies: sterol ester
"""
from rdkit import Chem

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    A sterol ester is a steroid ester obtained by formal condensation of the carboxy group of any carboxylic acid with the 3-hydroxy group of a sterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define steroid nucleus pattern (sterane backbone)
    steroid_nucleus_smiles = 'C1CCC2C(C1)CCC3C2CCC4C3CCCC4'  # Sterane skeleton
    steroid_nucleus = Chem.MolFromSmiles(steroid_nucleus_smiles)
    if steroid_nucleus is None:
        return False, "Error parsing steroid nucleus pattern"

    if not mol.HasSubstructMatch(steroid_nucleus):
        return False, "No steroid nucleus found"

    # Define ester functional group pattern
    ester_pattern = Chem.MolFromSmarts('C(=O)O')
    if ester_pattern is None:
        return False, "Error parsing ester pattern"

    # Find ester functional groups
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester group found"

    # Find steroid nucleus atom indices
    steroid_matches = mol.GetSubstructMatches(steroid_nucleus)
    if not steroid_matches:
        return False, "No steroid nucleus found"

    steroid_atom_indices = set()
    for match in steroid_matches:
        steroid_atom_indices.update(match)

    # Check if any ester oxygen is connected to the steroid nucleus
    ester_bonded_to_steroid = False
    for match in ester_matches:
        # Get the oxygen atom index in the ester group
        # The ester pattern is C(=O)O, so the third atom is oxygen
        oxygen_atom_idx = match[2]

        # Get the atom connected to the ester oxygen (excluding the carbon in the ester group)
        oxygen_atom = mol.GetAtomWithIdx(oxygen_atom_idx)
        neighbors = [nbr for nbr in oxygen_atom.GetNeighbors() if nbr.GetIdx() != match[1]]
        for neighbor in neighbors:
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx in steroid_atom_indices:
                ester_bonded_to_steroid = True
                break
        if ester_bonded_to_steroid:
            break

    if not ester_bonded_to_steroid:
        return False, "Ester group not connected to steroid nucleus"

    return True, "Contains steroid nucleus with esterified group connected"