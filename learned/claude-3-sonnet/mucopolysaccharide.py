"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is a mucopolysaccharide based on its SMILES string.
    Mucopolysaccharides are polysaccharides composed of alternating units from 
    uronic acids and glycosamines, commonly partially esterified with sulfuric acid.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a mucopolysaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More flexible patterns for sugar-like rings
    sugar_ring_pattern = Chem.MolFromSmarts("[C;R][C;R][C;R][C;R][C;R]")
    
    # Patterns for uronic acid variations (including lactones and modified forms)
    uronic_patterns = [
        Chem.MolFromSmarts("[C;R]C(=O)O"), # Simple uronic acid
        Chem.MolFromSmarts("[C;R]1[O;R]C(=O)"), # Lactone form
        Chem.MolFromSmarts("[C;R]C(=O)OC"), # Methyl ester
        Chem.MolFromSmarts("[C;R]C(=O)N"), # Uronic acid amide
    ]
    
    # Patterns for glycosamine variations
    glycosamine_patterns = [
        Chem.MolFromSmarts("[C;R][NH2,NH]"), # Simple amine
        Chem.MolFromSmarts("[C;R]NC(=O)C"), # N-acetyl
        Chem.MolFromSmarts("[C;R]N[C;R]"), # Ring nitrogen
        Chem.MolFromSmarts("[C;R]NC(=O)"), # Amide
    ]
    
    # Patterns for glycosidic and related linkages
    linkage_patterns = [
        Chem.MolFromSmarts("[C;R]O[C;R]"), # Classical glycosidic
        Chem.MolFromSmarts("[C;R]OC(=O)[C;R]"), # Ester linkage
        Chem.MolFromSmarts("[C;R]N[C;R]"), # N-glycosidic
    ]
    
    # Count matches
    sugar_rings = len(mol.GetSubstructMatches(sugar_ring_pattern))
    
    uronic_matches = sum(1 for pat in uronic_patterns if mol.HasSubstructMatch(pat))
    glycosamine_matches = sum(1 for pat in glycosamine_patterns if mol.HasSubstructMatch(pat))
    linkage_matches = sum(1 for pat in linkage_patterns if mol.HasSubstructMatch(pat))
    
    # Count rings
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    
    # Count key atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Scoring system
    score = 0
    if sugar_rings >= 1: score += 1
    if uronic_matches >= 1: score += 1
    if glycosamine_matches >= 1: score += 1
    if linkage_matches >= 1: score += 1
    if ring_count >= 2: score += 1
    if o_count >= 4: score += 1
    if n_count >= 1: score += 1
    if mol_wt >= 250: score += 1

    # Classification logic
    if score >= 6:
        reason = "Contains sugar-like rings with uronic acid and glycosamine components"
        if uronic_matches > 1 or glycosamine_matches > 1:
            reason += f" (found {uronic_matches} uronic and {glycosamine_matches} glycosamine units)"
        return True, reason
        
    # Specific rejection reasons
    if sugar_rings == 0:
        return False, "No sugar-like ring structures found"
    if uronic_matches == 0:
        return False, "No uronic acid components found"
    if glycosamine_matches == 0:
        return False, "No glycosamine components found"
    if linkage_matches == 0:
        return False, "No appropriate linkages found"
    
    return False, f"Insufficient structural features (score: {score}/8)"