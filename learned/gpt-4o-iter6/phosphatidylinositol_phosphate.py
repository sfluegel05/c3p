"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
from rdkit import Chem

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a member of the phosphatidylinositol phosphate class
    based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylinositol phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sn-glycerol backbone with two fatty acid esters (exclude inositol)
    glycerol_pattern = Chem.MolFromSmarts("C(COC(=O)[C@@H](CO)OC(=O)C)O")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with ester bonds found"
        
    # Look for phosphate group linked to glycerol backbone
    phosphate_glycerol_pattern = Chem.MolFromSmarts("COP([O-])(=O)[O-]")
    if not mol.HasSubstructMatch(phosphate_glycerol_pattern):
        return False, "No phosphate group linked to glycerol found"
        
    # Check for inositol ring with possible phosphates
    inositol_pattern = Chem.MolFromSmarts("C1(C(O)[CH]([CH]([CH]([CH](C1O)O)O)O)O)O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"

    # Count phosphate groups on inositol (at least one should be present)
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 1:
        return False, "No phosphate groups on inositol"

    # Successful classification as phosphatidylinositol phosphate
    return True, "Contains glycerol backbone, phosphates, and phosphorylated inositol ring"