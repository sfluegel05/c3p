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
    
    # Look for a more general glycerol backbone connected to fatty acids
    glycerol_pattern = Chem.MolFromSmarts("[C@@H](COC(=O)C)[O]C")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No flexible glycerol backbone with ester bonds found"

    # Look for phosphate group potentially linked to the glycerol or inositol
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) == 0:
        return False, "No phosphate group found"

    # Check for inositol ring with flexibility for phosphate attachments
    inositol_pattern = Chem.MolFromSmarts("C1(O)C(O)C(O)C(O)C(O)C1O")
    inositol_matches = mol.GetSubstructMatches(inositol_pattern)
    if not inositol_matches:
        return False, "No inositol ring found"

    # Count phosphate groups on inositol; allow for at least one phosphate attachment
    phosphate_inositol_pattern = Chem.MolFromSmarts("c1cc(OP([O-])=O)cc(OP([O-])=O)c1O")
    if not mol.HasSubstructMatch(phosphate_inositol_pattern):
        return False, "No phosphates on inositol, or unrecognizable placements"

    # Successful classification
    return True, "Contains flexible glycerol backbone, phosphates, and a phosphorylated inositol ring"