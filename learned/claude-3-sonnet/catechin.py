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
    Catechins are flavan-3-ols with various substitution patterns.
    
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

    # More flexible benzopyran core pattern that allows for substitutions
    # Matches the basic C6-C3-C6 skeleton with oxygen bridge
    core_pattern = Chem.MolFromSmarts("[#6]1-[#6]-[#6]-c2c(-[#6]1)c([#6,#8,#1])c([#6,#8,#1])c([#6,#8,#1])c2")
    if not mol.HasSubstructMatch(core_pattern):
        return False, "No flavonoid core structure found"

    # Check for the characteristic flavan-3-ol structure with more flexibility
    # Allows for various substitutions and both stereochemistries
    flavan3ol_pattern = Chem.MolFromSmarts("O1[C][C]([OH1])c2c(-[C]1)cccc2")
    if not mol.HasSubstructMatch(flavan3ol_pattern):
        return False, "No flavan-3-ol structure found"

    # Count hydroxyl groups (both aliphatic and aromatic)
    oh_pattern = Chem.MolFromSmarts("[OH1]")
    oh_matches = len(mol.GetSubstructMatches(oh_pattern))
    if oh_matches < 2:
        return False, f"Too few hydroxyl groups ({oh_matches}) for a catechin"

    # Look for characteristic substitution patterns
    features = []
    
    # Check for common A-ring substitution patterns (more flexible)
    a_ring_pattern = Chem.MolFromSmarts("O1[C][C]c2c(-[C]1)c(O)cc(O)c2")
    if mol.HasSubstructMatch(a_ring_pattern):
        features.append("5,7-dihydroxy pattern")

    # Check for B-ring patterns with more flexibility
    catechol_pattern = Chem.MolFromSmarts("c1c(O)c(O)cc([C])c1")
    pyrogallol_pattern = Chem.MolFromSmarts("c1c(O)c(O)c(O)cc1")
    
    if mol.HasSubstructMatch(pyrogallol_pattern):
        features.append("trihydroxyphenyl group")
    elif mol.HasSubstructMatch(catechol_pattern):
        features.append("dihydroxyphenyl group")

    # Check for common modifications
    ester_pattern = Chem.MolFromSmarts("OC(=O)")
    methoxy_pattern = Chem.MolFromSmarts("OC")
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)O")
    
    if mol.HasSubstructMatch(ester_pattern):
        features.append("ester derivative")
    if mol.HasSubstructMatch(methoxy_pattern):
        features.append("methoxy substitution")
    if mol.HasSubstructMatch(sulfate_pattern):
        features.append("sulfate group")

    # Basic structural checks
    ring_count = Chem.rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 2:
        return False, "Insufficient ring count for a catechin structure"

    mol_wt = Chem.Descriptors.ExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 1500:
        return False, f"Molecular weight ({mol_wt:.1f}) outside typical range for catechins"

    # Count aromatic rings
    arom_pattern = Chem.MolFromSmarts("a1aaaaa1")
    arom_rings = len(mol.GetSubstructMatches(arom_pattern))
    if arom_rings < 1:
        return False, "Missing required aromatic ring system"

    feature_str = ", ".join(features) if features else "basic"
    return True, f"Contains flavan-3-ol core with {oh_matches} hydroxyl groups. Features: {feature_str}"