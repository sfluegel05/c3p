"""
Classifies: CHEBI:35342 17alpha-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_17alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17alpha-hydroxy steroid based on its SMILES string.
    A 17alpha-hydroxy steroid has a steroid nucleus with a hydroxyl group at position 17 in the alpha configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for steroid nucleus (simplified pattern)
    steroid_pattern = Chem.MolFromSmarts("[C@H]1CC[C@H]2[C@@H]3CC[C@H]4C[C@H](O)CC[C@]4(C)[C@H]3CC[C@]12C")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid nucleus detected"

    # Find hydroxyl group(s)
    hydroxyls = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1:  # Assuming -OH
            hydroxyls.append(atom)

    if not hydroxyls:
        return False, "No hydroxyl groups"

    # Check if any hydroxyl is attached to C17 in alpha configuration
    # Assuming C17 is part of the D-ring (five-membered) in the steroid nucleus
    # Use SMARTS to find O attached to a carbon in the D-ring with alpha configuration
    # Alpha configuration might correspond to a specific chiral tag
    alpha_hydroxyl_pattern = Chem.MolFromSmarts("[C@](C)(C)(O)CC")
    alpha_matches = mol.GetSubstructMatches(alpha_hydroxyl_pattern)
    if not alpha_matches:
        return False, "No alpha-hydroxyl group detected"

    # Check if the alpha-hydroxyl is on the D-ring
    d_ring_pattern = Chem.MolFromSmarts("[C]1[C][C][C][C]1")  # Five-membered ring (D-ring)
    d_ring_matches = mol.GetSubstructMatches(d_ring_pattern)
    d_ring_atoms = set()
    for match in d_ring_matches:
        d_ring_atoms.update(match)

    for match in alpha_matches:
        carbon_idx = match[0]
        if carbon_idx in d_ring_atoms:
            return True, "17-alpha-hydroxy group detected on D-ring"

    return False, "No 17-alpha-hydroxy group found"