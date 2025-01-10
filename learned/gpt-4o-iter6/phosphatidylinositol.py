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

    # Confirm inositol ring with 5 hydroxyl groups
    inositol_oh_pattern = Chem.MolFromSmarts("C1([C@@H]([C@H]([C@H]([C@@H]([C@H]1O)O)O)O)O)O") 
    if not mol.HasSubstructMatch(inositol_oh_pattern):
        return False, "No inositol ring with 5 hydroxyl groups found"
    
    # Verify glycerol presence esterified to phosphate
    glycerol_pattern = Chem.MolFromSmarts("[O][CH2]C([O])([CH2][O][P](=O)([O])[O])")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone connected to phosphate"

    # Ensure phosphatidyl is connected correctly through ester bonds
    phosphatidyl_pattern = Chem.MolFromSmarts("C(OC=O)C(OC=O)C(OP(=O)(O)O)O")
    if not mol.HasSubstructMatch(phosphatidyl_pattern):
        return False, "No phosphatidyl ester linkage pattern found"

    return True, "Consistent with a phosphatidylinositol: inositol esterified with glycerophosphate connected to phosphatidyl groups"