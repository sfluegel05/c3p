"""
Classifies: CHEBI:28874 phosphatidylinositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phosphatidylinositol(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol based on its SMILES string.
    A phosphatidylinositol contains a myo-inositol ring, a glycerophosphate group, and two fatty acid chains.

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

    # Look for myo-inositol ring pattern (six-membered carbon ring with five or six hydroxy groups)
    inositol_pattern = Chem.MolFromSmarts("C1(C(C(C(C(C1O)O)O)O)O)O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No myo-inositol ring structure found"

    # Look for phosphodiester bond linking the inositol ring
    phosphoester_pattern = Chem.MolFromSmarts("COP(=O)(O)O")
    if not mol.HasSubstructMatch(phosphoester_pattern):
        return False, "No phosphate group attached to inositol ring"

    # Check for glycerol backbone with two esterified fatty acid chains
    glycerol_fatty_acid_pattern = Chem.MolFromSmarts("OCC(=O)CC(=O)O")
    glycerol_matches = mol.GetSubstructMatches(glycerol_fatty_acid_pattern)
    if len(glycerol_matches) < 1:
        return False, "No glycerophospholipid structure (glycerol with two fatty acid chains) found"

    return True, "Contains myo-inositol ring with glycerophosphate and two esterified fatty acid chains"