"""
Classifies: CHEBI:50699 oligosaccharide
"""
"""
Classifies: oligosaccharide
A compound in which monosaccharide units are joined by glycosidic linkages.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_oligosaccharide(smiles: str):
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an oligosaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for sugar rings (6-membered oxane or 5-membered oxolane)
    oxane_pattern = Chem.MolFromSmarts("[C]1[C][C][C][C]O1")  # 6-membered sugar ring
    oxolane_pattern = Chem.MolFromSmarts("[C]1[C][C][C]O1")   # 5-membered sugar ring
    
    oxane_matches = len(mol.GetSubstructMatches(oxane_pattern))
    oxolane_matches = len(mol.GetSubstructMatches(oxolane_pattern))
    total_rings = oxane_matches + oxolane_matches
    
    if total_rings < 2:
        return False, f"Found only {total_rings} sugar rings, need at least 2"
    
    if total_rings > 20:
        return False, f"Found {total_rings} sugar rings, likely a polysaccharide"
        
    # Look for glycosidic linkages (C-O-C between rings)
    glycosidic_pattern = Chem.MolFromSmarts("[C][O][C]")
    glycosidic_matches = len(mol.GetSubstructMatches(glycosidic_pattern))
    
    if glycosidic_matches < 1:
        return False, "No glycosidic linkages found"
        
    # Check for hydroxyl groups (characteristic of sugars)
    hydroxyl_pattern = Chem.MolFromSmarts("[O][H]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    
    if hydroxyl_matches < 3:
        return False, "Too few hydroxyl groups for an oligosaccharide"
        
    # Count oxygen atoms (should be abundant in sugars)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 5:
        return False, "Too few oxygen atoms for an oligosaccharide"
        
    # Calculate molecular weight (should be within reasonable range for oligosaccharides)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, "Molecular weight too low for oligosaccharide"
    if mol_wt > 5000:
        return False, "Molecular weight too high, likely a polysaccharide"
        
    # Check ring connectivity
    ring_info = mol.GetRingInfo()
    if not ring_info.NumRings():
        return False, "No rings found"
        
    # Additional check for characteristic sugar features
    sugar_pattern = Chem.MolFromSmarts("[CH1,CH2][OH]")  # CHOH or CH2OH groups
    sugar_matches = len(mol.GetSubstructMatches(sugar_pattern))
    if sugar_matches < 2:
        return False, "Missing characteristic sugar hydroxyl groups"

    return True, f"Contains {total_rings} sugar rings connected by glycosidic linkages with appropriate hydroxyl groups"