"""
Classifies: CHEBI:18179 phosphoinositide
"""
from rdkit import Chem

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    A phosphoinositide is a phosphatidylinositol with one or more phosphate groups on the inositol ring.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        (bool, str): Tuple where first element indicates if the molecule is a phosphoinositide,
                      and the second element provides the reason for the classification.
    """
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the pattern for the myo-inositol ring (more generic)
    inositol_ring_pattern = Chem.MolFromSmarts("C1(C(C(C(C(C1O)O)O)O)O)O")
    if not mol.HasSubstructMatch(inositol_ring_pattern):
        return False, "No full inositol ring found"

    # Define the pattern for phosphate groups attached to an alcohol group
    phosphate_group_pattern = Chem.MolFromSmarts("OP(=O)(O)[O-]")
    
    # Check for phosphates attached to inositol ring
    phosphate_matches = mol.GetSubstructMatches(phosphate_group_pattern)
    phosphate_count = len(phosphate_matches)
    
    # At least one phosphate group must be present
    if phosphate_count < 1:
        return False, f"Phosphate groups insufficient: found {phosphate_count}, need at least 1"

    # Define a pattern to validate phosphatidyl backbone (basic pattern for general glycerol-phosphate linkage)
    glycerol_linkage_pattern = Chem.MolFromSmarts("C(O)C(COC(=O)C)OP(=O)(O)O")
    if not mol.HasSubstructMatch(glycerol_linkage_pattern):
        return False, "Phosphatidyl linkage not found"

    return True, "Molecule has a phosphatidylinositol backbone with phosphate attachments on the inositol ring"