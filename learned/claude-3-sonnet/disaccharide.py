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

    # Count rings
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    if ring_count > 3:  # Allow up to 3 rings for special cases
        return False, f"Too many rings ({ring_count}), disaccharides typically have 2"
    if ring_count < 1:
        return False, f"Too few rings ({ring_count}), disaccharides need at least 1"

    # Look for sugar ring patterns more specifically
    pyranose = Chem.MolFromSmarts("[C]1[C][C][C]([C][O]1)")
    furanose = Chem.MolFromSmarts("[C]1[C][C]([C][O]1)")
    kdo = Chem.MolFromSmarts("CC(=O)[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO")
    
    sugar_ring_matches = (len(mol.GetSubstructMatches(pyranose)) + 
                         len(mol.GetSubstructMatches(furanose)) +
                         len(mol.GetSubstructMatches(kdo)))
    
    if sugar_ring_matches < 1:
        return False, "No sugar ring patterns found"

    # Check for glycosidic bond pattern with more specific context
    glycosidic_pattern = Chem.MolFromSmarts("[C][OX2][C;R]")  # Require one carbon to be in ring
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    if len(glycosidic_matches) < 1:
        return False, "No glycosidic bond found"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Expanded ranges to accommodate modified sugars and special cases
    if c_count < 10 or c_count > 18:  # Allow for Kdo sugars
        return False, f"Carbon count ({c_count}) outside range for disaccharides (10-18)"
    
    if o_count < 8 or o_count > 14:  # Allow for uronic acids and modifications
        return False, f"Oxygen count ({o_count}) outside range for disaccharides (8-14)"

    # Look for hydroxyl groups with more context
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OX2H1]")  # Require sp3 carbon
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 4:  # Reduced minimum to account for modifications
        return False, f"Too few hydroxyl groups ({len(hydroxyl_matches)}), need at least 4"

    # Check for characteristic sugar carbon patterns
    sugar_carbon_pattern = Chem.MolFromSmarts("[C]([O])([C])[C]([O])")
    sugar_carbon_matches = mol.GetSubstructMatches(sugar_carbon_pattern)
    if len(sugar_carbon_matches) < 4:
        return False, "Insufficient sugar carbon pattern matches"

    # Check molecular weight with expanded range
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 800:  # Expanded range for modified disaccharides
        return False, f"Molecular weight ({mol_wt}) outside range for disaccharides (250-800)"

    # Additional check for non-sugar atoms
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    if n_count > 2 or s_count > 1:  # Allow limited N for amino sugars
        return False, "Too many non-sugar atoms"

    # Check for characteristic disaccharide connectivity
    sugar_connection = Chem.MolFromSmarts("[C;R]([O])([C])[O][C;R]([O])([C])")
    if not mol.HasSubstructMatch(sugar_connection):
        return False, "Missing characteristic sugar-sugar connection"

    return True, "Contains two sugar units connected by glycosidic bond with appropriate number of hydroxyl groups"