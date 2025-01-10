"""
Classifies: CHEBI:23053 catechin
"""
"""
Classifies: CHEBI:23053 catechin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.
    Catechins are flavan-3-ols with a characteristic structure and substitution pattern.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a catechin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic flavan-3-ol core (more flexible pattern)
    # Matches the benzopyran core with 3-OH, allowing for various substitutions
    core_pattern = Chem.MolFromSmarts("[OX2][CH1]1[CH2][c]2c([CH1]1)cccc2")
    
    # Alternative core pattern to catch more complex derivatives
    core_pattern2 = Chem.MolFromSmarts("O1[CH1][CH2]c2c([CH1]1)c([OH1,O])cc([OH1,O])c2")
    
    if not (mol.HasSubstructMatch(core_pattern) or mol.HasSubstructMatch(core_pattern2)):
        return False, "Missing required flavan-3-ol core structure"

    # Count rings to ensure basic structure
    ring_count = Chem.rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 2:
        return False, "Insufficient ring count for catechin structure"

    # Check for characteristic A-ring substitution (5,7-dihydroxy pattern)
    a_ring = Chem.MolFromSmarts("O1[CH1][CH2]c2c([CH1]1)c([OH1,O])cc([OH1,O])c2")
    
    # Check for B-ring patterns (allowing for various hydroxylation patterns)
    b_ring_dihydroxy = Chem.MolFromSmarts("c1c([OH1,O])c([OH1,O])ccc1")
    b_ring_trihydroxy = Chem.MolFromSmarts("c1c([OH1,O])c([OH1,O])c([OH1,O])cc1")
    
    features = []
    
    if mol.HasSubstructMatch(a_ring):
        features.append("5,7-dihydroxy pattern")
    
    if mol.HasSubstructMatch(b_ring_trihydroxy):
        features.append("trihydroxyphenyl B-ring")
    elif mol.HasSubstructMatch(b_ring_dihydroxy):
        features.append("dihydroxyphenyl B-ring")
        
    # Check for common modifications
    galloyl_pattern = Chem.MolFromSmarts("C(=O)c1c([OH1])c([OH1])c([OH1])cc1")
    methoxy_pattern = Chem.MolFromSmarts("OC")
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)[OH1]")
    
    if mol.HasSubstructMatch(galloyl_pattern):
        features.append("galloyl group")
    if mol.HasSubstructMatch(methoxy_pattern):
        features.append("methoxy substitution")
    if mol.HasSubstructMatch(sulfate_pattern):
        features.append("sulfate group")

    # Count hydroxyl groups
    oh_pattern = Chem.MolFromSmarts("[OH1]")
    oh_count = len(mol.GetSubstructMatches(oh_pattern))
    
    # Basic structural requirements
    if oh_count < 2:
        return False, f"Insufficient hydroxyl groups ({oh_count}) for catechin"
    
    mol_wt = Chem.Descriptors.ExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 2000:
        return False, f"Molecular weight ({mol_wt:.1f}) outside typical range for catechins"

    # Additional check for characteristic oxygen pattern
    o_pattern = Chem.MolFromSmarts("[OX2]")
    o_count = len(mol.GetSubstructMatches(o_pattern))
    if o_count < 3:
        return False, "Insufficient oxygen atoms for catechin structure"

    # If we have the core and at least some characteristic features, classify as catechin
    if features:
        feature_str = ", ".join(features)
        return True, f"Catechin derivative with {feature_str}. Contains {oh_count} hydroxyl groups"
    
    # If we have the core but no characteristic features, be more cautious
    return False, "Has basic structure but lacks characteristic catechin substitution patterns"