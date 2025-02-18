"""
Classifies: CHEBI:18179 phosphoinositide
"""
"""
Classifies: CHEBI:18179 phosphoinositide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide (phosphorylated phosphatidylinositol).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if phosphoinositide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Step 1: Verify phosphatidylinositol core structure
    # Glycerol backbone with two ester groups and a phosphate linked to inositol
    # Inositol is a 6-membered carbocycle with multiple oxygen substituents
    # SMARTS pattern allows variable substitution on inositol (OH or phosphate groups)
    pi_core_pattern = Chem.MolFromSmarts(
        "[CH2]-[CH](-[OX2]-C(=O)[!O])-[CH2]-[OX2]-P(=[OX1])(-[OX2])-[OX2][C@H]1[C@H]([OX2,OX1,C(=O),P])[C@H]([OX2,OX1,C(=O),P])[C@H]([OX2,OX1,C(=O),P])[C@H]([OX2,OX1,C(=O),P])[C@H]1[OX2,OX1,C(=O),P]"
    )
    if not mol.HasSubstructMatch(pi_core_pattern):
        return False, "Not a phosphatidylinositol"
    
    # Step 2: Check for at least one additional phosphate group on inositol
    # Look for O-P(=O)(O) groups attached to the inositol ring (excluding the connecting phosphate)
    # The connecting phosphate is part of the glycerol-phospho-inositol linkage
    phosphate_pattern = Chem.MolFromSmarts("[C][OX2]-P(=[OX1])(-[OX2])-[OX2]")
    matches = mol.GetSubstructMatches(phosphate_pattern)
    
    # The core structure includes one phosphate (glycerol-inositol link), so need at least two in total
    if len(matches) < 2:
        return False, f"Only {len(matches)} phosphate groups found (needs ≥2)"
    
    return True, "Phosphatidylinositol with additional phosphorylation"