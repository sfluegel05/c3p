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
    
    # Look for sugar rings with proper substitution patterns
    # Pyranose (6-membered) ring pattern with typical sugar substitutions
    pyranose_pattern = Chem.MolFromSmarts("[C]1[C]([OH0,OH1])[C]([OH0,OH1])[C]([OH0,OH1])[C]([OH0,OH1])O1")
    # Furanose (5-membered) ring pattern with typical sugar substitutions
    furanose_pattern = Chem.MolFromSmarts("[C]1[C]([OH0,OH1])[C]([OH0,OH1])[C]([OH0,OH1])O1")
    
    pyranose_matches = len(mol.GetSubstructMatches(pyranose_pattern))
    furanose_matches = len(mol.GetSubstructMatches(furanose_pattern))
    total_rings = pyranose_matches + furanose_matches
    
    if total_rings < 2:
        return False, f"Found only {total_rings} sugar rings, need at least 2"
    
    if total_rings > 20:
        return False, f"Found {total_rings} sugar rings, likely a polysaccharide"
        
    # Look for glycosidic linkages between rings
    # More specific pattern that looks for C-O-C between ring carbons
    glycosidic_pattern = Chem.MolFromSmarts("[C;R][O][C;R]")
    glycosidic_matches = len(mol.GetSubstructMatches(glycosidic_pattern))
    
    if glycosidic_matches < 1:
        return False, "No glycosidic linkages found"
        
    # Check for hydroxyl and hydroxyl-derived groups (including modified ones)
    hydroxyl_pattern = Chem.MolFromSmarts("[$([OH1]),$([OH0][C]=[O]),$(O[CH3]),$(O[C]=[O])]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    
    if hydroxyl_matches < 4:
        return False, "Too few hydroxyl or hydroxyl-derived groups"
        
    # Count oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 6:
        return False, "Too few oxygen atoms for an oligosaccharide"
        
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:  # Lowered threshold
        return False, "Molecular weight too low for oligosaccharide"
    if mol_wt > 10000:  # Increased threshold
        return False, "Molecular weight too high, likely a polysaccharide"
        
    # Check for characteristic sugar features
    # Look for CH2OH groups (common in sugars)
    ch2oh_pattern = Chem.MolFromSmarts("[CH2][OH1]")
    # Look for CHOH groups
    choh_pattern = Chem.MolFromSmarts("[CH1][OH1]")
    # Look for modified versions (acetylated, methylated, etc)
    modified_pattern = Chem.MolFromSmarts("[CH1,CH2][O][C]")
    
    sugar_features = (
        len(mol.GetSubstructMatches(ch2oh_pattern)) +
        len(mol.GetSubstructMatches(choh_pattern)) +
        len(mol.GetSubstructMatches(modified_pattern))
    )
    
    if sugar_features < 3:
        return False, "Missing characteristic sugar structural features"

    # Check for reasonable C:O ratio for sugars
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count > 0 and o_count/c_count < 0.3:
        return False, "Carbon to oxygen ratio not consistent with oligosaccharides"

    return True, f"Contains {total_rings} sugar rings connected by glycosidic linkages with appropriate sugar characteristics"