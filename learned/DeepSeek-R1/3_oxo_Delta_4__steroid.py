"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
"""
Classifies: 3-oxo-Delta(4) steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES string.
    A 3-oxo-Delta(4) steroid has a ketone at position 3 and a conjugated double bond at C4-C5 in a steroid nucleus.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule matches criteria, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check steroid nucleus using SMARTS pattern for tetracyclic system
    steroid_smarts = Chem.MolFromSmarts("[C@]12[C@]3[C@]([C@H]([C@@H]1CC2)C)(CC[C@H]4[C@@]3(CC[C@H](C4)O)C)C")
    if not mol.HasSubstructMatch(steroid_smarts):
        return False, "No steroid nucleus detected"

    # Locate 3-oxo group in ring A (should be in 6-membered ring)
    ring_a_ketone = Chem.MolFromSmarts("[C;R1](=O)[C;R1][C;R1][C;R1][C;R1][C;R1]")
    ketone_matches = mol.GetSubstructMatches(ring_a_ketone)
    if not ketone_matches:
        return False, "3-oxo group not found in ring A"

    # Verify conjugated double bond in positions 4-5 (adjacent to ketone in same ring)
    conjugated_system = Chem.MolFromSmarts("[C;R1]=[C;R1]-[C;R1](=O)")
    conjugated_matches = mol.GetSubstructMatches(conjugated_system)
    if not conjugated_matches:
        return False, "No conjugated Delta(4) system adjacent to 3-oxo group"

    # Molecular weight sanity check (typical steroids >280 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 280:
        return False, f"Molecular weight ({mol_wt:.1f} Da) too low for steroid"

    return True, "Contains 3-oxo group and conjugated Delta(4) system in steroid nucleus"