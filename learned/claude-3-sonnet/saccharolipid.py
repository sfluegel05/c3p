"""
Classifies: CHEBI:166828 saccharolipid
"""
"""
Classifies: saccharolipids - lipids containing carbohydrate moieties
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    Saccharolipids are lipids that contain carbohydrate moieties with fatty acids 
    directly attached to the sugar units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saccharolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sugar rings (pyranose/furanose) with specific hydroxylation pattern
    sugar_pattern = Chem.MolFromSmarts("[CR1]1[CR1]([OH1,OR])[CR1]([OH1,OR])[CR1]([OH1,OR])[CR1][OR1]1")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No appropriately hydroxylated sugar rings found"

    # Look for fatty acid chains
    fatty_acid_pattern = Chem.MolFromSmarts("[CH2][CH2][CH2][CH2][CH2][CH2]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    
    # Look for acylation (fatty acid attachment to sugar)
    acyl_sugar_pattern = Chem.MolFromSmarts("[CR1]1[CR1][CR1][CR1][CR1][OR1]1-[CH2]-O-C(=O)-[CH2][CH2][CH2]")
    direct_acylation = mol.HasSubstructMatch(acyl_sugar_pattern)
    
    if not (fatty_acid_matches and direct_acylation):
        return False, "No direct fatty acid attachment to sugar found"

    # Count hydroxyl groups on sugar rings
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    
    # Look for characteristic linkages
    ester_pattern = Chem.MolFromSmarts("[#6]-C(=O)-O-[#6]")
    has_ester = mol.HasSubstructMatch(ester_pattern)
    
    if not has_ester:
        return False, "No ester linkages found"
        
    # Optional: Check for phosphate groups (common in many saccharolipids)
    phosphate_pattern = Chem.MolFromSmarts("[P](=O)([O-,OH])([O-,OH])")
    has_phosphate = mol.HasSubstructMatch(phosphate_pattern)
    
    # Count key elements
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    num_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if num_oxygens < 4:
        return False, "Insufficient oxygen atoms for saccharolipid"

    # Check for characteristic saccharolipid features
    if len(sugar_matches) >= 1 and len(fatty_acid_matches) >= 1 and direct_acylation:
        reason = f"Contains {len(sugar_matches)} sugar ring(s) with direct fatty acid attachment, "
        reason += f"{len(hydroxyl_matches)} hydroxyl groups, "
        if has_phosphate:
            reason += "phosphate groups, "
        reason += "and appropriate ester linkages"
        return True, reason
        
    return False, "Does not meet structural requirements for saccharolipid classification"