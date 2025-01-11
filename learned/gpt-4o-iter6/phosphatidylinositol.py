"""
Classifies: CHEBI:28874 phosphatidylinositol
"""
from rdkit import Chem

def is_phosphatidylinositol(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol based on its SMILES string.
    A phosphatidylinositol is a glycerophosphoinositol having a phosphatidyl group esterified to one of the hydroxy groups of inositol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylinositol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for inositol ring (a six-membered carbon ring with 5 or 6 hydroxyl groups)
    inositol_pattern = Chem.MolFromSmarts("C1(C(C(C(C(C1O)O)O)O)O)O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"

    # Check for phosphatidyl group (a phosphate group with ester linkages)
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(=O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphatidyl group found"
    
    # Check for glycerol backbone as part of the phosphatidyl moiety
    glycerol_pattern = Chem.MolFromSmarts("C(C(CO)O)O")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Confirm one of the inositol hydroxyl groups is esterified with the phosphatidyl group
    ester_linkage_pattern = Chem.MolFromSmarts("[C@H](O[*:1])[#6]")
    if not mol.HasSubstructMatch(ester_linkage_pattern):
        return False, "No ester linkage between phosphatidyl group and inositol found"
    
    return True, "Contains glycerol backbone with phosphatidyl group esterified to inositol"