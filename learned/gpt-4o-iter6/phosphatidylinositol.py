"""
Classifies: CHEBI:28874 phosphatidylinositol
"""
from rdkit import Chem

def is_phosphatidylinositol(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol based on its SMILES string.
    A phosphatidylinositol has a phosphatidyl group esterified to one of the hydroxy groups of inositol.

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

    # Check of inositol scaffold is present with 5 hydroxyl groups
    inositol_oh_pattern = Chem.MolFromSmarts("C1([C@@H]([C@H]([C@H]([C@@H]([C@H]1O)O)O)O)O)O")
    if not mol.HasSubstructMatch(inositol_oh_pattern):
        return False, "No inositol ring with 5 hydroxyl groups found"

    # Look for the phosphatidyl group ester linkage
    phosphatidyl_pattern = Chem.MolFromSmarts("C(OC(=O)C)OP(=O)(O)O")
    if not mol.HasSubstructMatch(phosphatidyl_pattern):
        return False, "No phosphatidyl ester linkage found"
        
    # Check for full glycerol presence in phosphatidyl group
    glycerol_pattern = Chem.MolFromSmarts("C(C(CO)O)O")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found in association with phosphate"

    # Ensure linkage of glycerol to inositol
    esterified_inositol_pattern = Chem.MolFromSmarts("OC[C@H](O[CH2])OP")
    if not mol.HasSubstructMatch(esterified_inositol_pattern):
        return False, "No ester linkage correctly connecting glycerol phosphate to inositol"

    return True, "Contains glycerol backbone with phosphatidyl group esterified to inositol"