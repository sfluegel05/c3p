"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
"""
Classifies: 3beta-hydroxy-Delta(5)-steroid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy-Delta(5)-steroid based on SMILES.
    Criteria:
    - Contains steroid nucleus (four fused rings: 3x six-membered, 1x five-membered)
    - 3beta-hydroxyl group (OH on carbon 3 in beta configuration)
    - Double bond between positions 5 and 6 in the steroid structure
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES"

    # Steroid nucleus check (crude approximation)
    steroid_core = Chem.MolFromSmarts("[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid nucleus detected"

    # Check for 3beta-hydroxyl group
    # Carbon 3 is in ring A, beta configuration (same face as C18/C19 methyl groups)
    hydroxyl_pattern = Chem.MolFromSmarts("[C@H]1[C@@H](O)CC[C@H]2[C@@H]1CC[C@@H]1[C@@H]2CCCC1")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "3beta-hydroxyl group not found"

    # Check for 5-6 double bond in the B ring
    double_bond_pattern = Chem.MolFromSmarts("[C@@H]1C=C[C@@H]2")
    matches = mol.GetSubstructMatches(double_bond_pattern)
    if not matches:
        return False, "No 5-6 double bond found"

    return True, "3beta-hydroxy-steroid with Delta(5) double bond"