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

    # More flexible uronic acid pattern (sugar ring with carboxylic acid)
    uronic_acid_pattern = Chem.MolFromSmarts("[C;R][C;R][C;R][C;R][C;R]C(=O)O")
    
    # More flexible glycosamine pattern (sugar ring with amine)
    glycosamine_pattern = Chem.MolFromSmarts("[C;R][C;R][C;R][C;R][C;R][NH2,NH]")
    
    # Glycosidic linkage pattern
    glycosidic_pattern = Chem.MolFromSmarts("[C;R]O[C;R]")
    
    # Sulfate ester pattern
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)O")
    
    # Check for minimum structural requirements
    uronic_matches = len(mol.GetSubstructMatches(uronic_acid_pattern))
    glycosamine_matches = len(mol.GetSubstructMatches(glycosamine_pattern))
    glycosidic_matches = len(mol.GetSubstructMatches(glycosidic_pattern))
    has_sulfate = mol.HasSubstructMatch(sulfate_pattern)
    
    # Count rings (should have multiple sugar rings)
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    
    # Basic requirements checks
    if ring_count < 2:
        return False, "Insufficient ring structures for a polysaccharide"
    
    if uronic_matches < 1:
        return False, "No uronic acid units found"
        
    if glycosamine_matches < 1:
        return False, "No glycosamine units found"
        
    if glycosidic_matches < 1:
        return False, "No glycosidic linkages found"
    
    # Count key atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if o_count < 6:
        return False, "Too few oxygen atoms for a polysaccharide"
        
    if n_count < 1:
        return False, "No nitrogen atoms found (required for glycosamines)"
    
    # Check molecular weight (should be substantial for a polysaccharide)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:  # Adjusted threshold
        return False, "Molecular weight too low for a mucopolysaccharide"
    
    # Success conditions
    success_msg = "Contains uronic acid and glycosamine units with glycosidic linkages"
    if has_sulfate:
        success_msg += " and sulfate groups"
    
    # Additional check for reasonable ratio of components
    if uronic_matches >= 1 and glycosamine_matches >= 1 and glycosidic_matches >= 1:
        return True, success_msg
    
    return False, "Does not have proper balance of structural components"