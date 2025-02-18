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

    # Check for 3-oxo group (carbonyl in a ring)
    ketone_pattern = Chem.MolFromSmarts('[C;R]=[O]')
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "No 3-oxo group detected"

    # Check for conjugated double bond at Delta(4) position (C=C adjacent to ketone in same ring)
    conjugated_pattern = Chem.MolFromSmarts('[C;R](=[O])-[C;R]=[C;R]')
    conjugated_matches = mol.GetSubstructMatches(conjugated_pattern)
    if not conjugated_matches:
        return False, "No conjugated double bond at alpha,beta position"

    # Verify steroid nucleus characteristics (4+ fused rings)
    ring_info = mol.GetRingInfo()
    if len(ring_info.AtomRings()) < 4:
        return False, "Insufficient rings for steroid structure"

    # Basic molecular weight check (steroids typically >250 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, f"Molecular weight ({mol_wt:.1f} Da) too low for steroid"

    return True, "Contains 3-oxo group and conjugated Delta(4) system in steroid nucleus"