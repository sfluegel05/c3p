"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
"""
Classifies: CHEBI:28875 myo-inositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    A myo-inositol phosphate is an inositol phosphate with a myo-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a myo-inositol phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for 6-membered ring (inositol core)
    ring_info = mol.GetRingInfo()
    if not any(len(ring) == 6 for ring in ring_info.AtomRings()):
        return False, "No 6-membered ring found"

    # Check for myo-inositol backbone pattern (6 carbons with hydroxyl/phosphate groups)
    # More flexible pattern that doesn't enforce specific stereochemistry
    myo_inositol_pattern = Chem.MolFromSmarts("[C]1([C]([OH,OP])([C]([OH,OP])([C]([OH,OP])([C]([OH,OP])([C]1[OH,OP])))))")
    if not mol.HasSubstructMatch(myo_inositol_pattern):
        return False, "No myo-inositol backbone pattern found"

    # Look for phosphate groups (more comprehensive pattern)
    phosphate_pattern = Chem.MolFromSmarts("[OX2][PX4](=[OX1])([OX2H,OX2-])[OX2H,OX2-]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 1:
        return False, "No phosphate groups found"

    # Check molecular formula (should have C6H(6-12)O(6-12)P(1-6))
    formula = rdMolDescriptors.CalcMolFormula(mol)
    c_count = formula.count('C')
    h_count = formula.count('H')
    o_count = formula.count('O')
    p_count = formula.count('P')
    
    if c_count != 6:
        return False, f"Wrong number of carbons: {c_count}"
    if h_count < 6 or h_count > 12:
        return False, f"Wrong number of hydrogens: {h_count}"
    if o_count < 6 or o_count > 12:
        return False, f"Wrong number of oxygens: {o_count}"
    if p_count < 1 or p_count > 6:
        return False, f"Wrong number of phosphates: {p_count}"

    return True, "Contains myo-inositol backbone with at least one phosphate group"