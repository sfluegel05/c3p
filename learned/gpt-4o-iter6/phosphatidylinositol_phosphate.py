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

    # Look for sn-glycerol backbone generally (flexibly considering positions of esters)
    glycerol_pattern = Chem.MolFromSmarts("[C@@H](COC(=O)[C,H])C(O)O(C=O[C,H])")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with ester bonds found"
        
    # Look for phosphate group linked to glycerol backbone
    phosphate_glycerol_pattern = Chem.MolFromSmarts("COP(=O)(O)")
    if not mol.HasSubstructMatch(phosphate_glycerol_pattern):
        return False, "No phosphate group linked to glycerol found"
        
    # Check for inositol ring with possible phosphates (pattern allows for flexibility)
    inositol_pattern = Chem.MolFromSmarts("C1(C[CH]([CH]([CH](C1O)O)O)O)O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found with potential phosphates"

    # Count phosphate groups on inositol (at least one should be present, considering flexibility)
    phosphate_pattern = Chem.MolFromSmarts("OP([O-])=O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 1:
        return False, "No phosphate groups on inositol"

    # Successful classification as phosphatidylinositol phosphate
    return True, "Contains glycerol backbone, phosphates, and a phosphorylated inositol ring"