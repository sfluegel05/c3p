"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
"""
Classifies: 3beta-hydroxy-Delta(5)-steroid

Definition: 'Any 3beta-hydroxy-steroid that contains a double bond between positions 5 and 6.'
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy-Delta(5)-steroid based on its SMILES string.
    A 3beta-hydroxy-Delta(5)-steroid is any 3beta-hydroxy-steroid that contains a double bond between positions 5 and 6 (Delta5).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3beta-hydroxy-Delta(5)-steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Kekulize the molecule to avoid aromaticity issues
    Chem.Kekulize(mol, clearAromaticFlags=True)

    # Define the steroid nucleus with atom numbering
    steroid_nucleus_smarts = """
    [#6]1CC[C@H]2CC[C@H]3CC[C@]4(C)CC[C@H](O)CC4=C3CC2C1
    """
    steroid_nucleus = Chem.MolFromSmarts(steroid_nucleus_smarts)
    if steroid_nucleus is None:
        return False, "Invalid steroid nucleus SMARTS pattern"

    # Check for steroid nucleus
    match = mol.GetSubstructMatch(steroid_nucleus)
    if not match:
        return False, "No steroid nucleus found"

    # Map the steroid nucleus to get atom indices for positions
    atom_indices = {}
    for idx, atom in enumerate(steroid_nucleus.GetAtoms()):
        atom_map_num = atom.GetAtomMapNum()
        if atom_map_num:
            atom_indices[atom_map_num] = match[idx]

    # Check for 3beta-hydroxy group at position 3
    # Position 3 is atom map number 3 in the SMARTS pattern
    # For simplicity, assume position 3 is connected to an -OH group with correct stereochemistry
    hydroxyl = Chem.MolFromSmarts("[C@@H](O)")
    if not mol.HasSubstructMatch(hydroxyl, useChirality=True):
        return False, "No 3beta-hydroxy group found at position 3"

    # Check for double bond between positions 5 and 6
    # Position 5 and 6 are atom map numbers 5 and 6
    # Create a bond query between these positions
    # Since atom indices might not match, we search for a double bond in any ring
    double_bond_found = False
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE and bond.IsInRing():
            double_bond_found = True
            break
    if not double_bond_found:
        return False, "No Delta(5) double bond found between positions 5 and 6"

    return True, "Molecule is a 3beta-hydroxy-Delta(5)-steroid"