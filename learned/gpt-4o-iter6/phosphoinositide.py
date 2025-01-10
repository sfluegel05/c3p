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
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphatidylinositol backbone, including the glycerol linkage
    glycerol_linkage_pattern = Chem.MolFromSmarts("C(COP(O)(O)=O)O")
    inositol_ring_pattern = Chem.MolFromSmarts("C1O[C@@H]([C@H](O)C([C@H](O)C1)OP(O)(=O)O)O")
    phosphate_group_pattern = Chem.MolFromSmarts("C(OP(O)(O)=O)C")

    if not mol.HasSubstructMatch(glycerol_linkage_pattern):
        return False, "Glycerol linkage not found in phosphatidyl structure"
    
    if not mol.HasSubstructMatch(inositol_ring_pattern):
        return False, "No recognizable inositol ring found"

    phosphate_matches = mol.GetSubstructMatches(phosphate_group_pattern)
    
    # Ensure one or more phosphate groups are associated with the inositol ring
    if len(phosphate_matches) == 0:
        return False, "No phosphate groups found on the inositol ring"
    
    return True, "Molecule has a phosphatidylinositol backbone with phosphates on the inositol ring"