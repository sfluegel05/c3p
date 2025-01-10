"""
Classifies: CHEBI:36233 disaccharide
"""
"""
Classifies: CHEBI:18667 disaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide consists of two monosaccharides joined by a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a disaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count rings - should have exactly 2 for most disaccharides
    # (some rare cases might have one ring opened)
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    if ring_count > 2:
        return False, f"Too many rings ({ring_count}), disaccharides typically have 2"
    if ring_count < 1:
        return False, f"Too few rings ({ring_count}), disaccharides typically have 2"

    # Check for glycosidic bond pattern (C-O-C between sugars)
    glycosidic_pattern = Chem.MolFromSmarts("[C][OX2][C]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    if len(glycosidic_matches) < 1:
        return False, "No glycosidic bond found"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Most common disaccharides have 12 carbons (6+6), but some have 11 (5+6) or 10 (5+5)
    if c_count < 10 or c_count > 13:
        return False, f"Carbon count ({c_count}) outside typical range for disaccharides (10-13)"
    
    # Typical O count for disaccharides (including ring O and OH groups)
    if o_count < 9 or o_count > 11:
        return False, f"Oxygen count ({o_count}) outside typical range for disaccharides (9-11)"

    # Look for hydroxyl group pattern typical in sugars
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 6:
        return False, f"Too few hydroxyl groups ({len(hydroxyl_matches)}), need at least 6"

    # Check for sugar ring pattern (pyranose or furanose)
    pyranose = Chem.MolFromSmarts("[C]1[C][C][C]([C][O]1)")
    furanose = Chem.MolFromSmarts("[C]1[C][C]([C][O]1)")
    
    sugar_ring_count = len(mol.GetSubstructMatches(pyranose)) + len(mol.GetSubstructMatches(furanose))
    if sugar_ring_count < 1:
        return False, "No sugar ring patterns found"

    # Additional check for carbohydrate-like structure
    # Most carbons should have a hydroxyl group attached
    carbons_with_oh = sum(1 for match in mol.GetSubstructMatches(Chem.MolFromSmarts("[C][OH]")))
    if carbons_with_oh < 6:
        return False, f"Too few carbons with hydroxyl groups ({carbons_with_oh}), need at least 6"

    # Check molecular weight - should be in typical range for disaccharides
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 280 or mol_wt > 400:
        return False, f"Molecular weight ({mol_wt}) outside typical range for disaccharides (280-400)"

    return True, "Contains two sugar units connected by glycosidic bond with appropriate number of hydroxyl groups"