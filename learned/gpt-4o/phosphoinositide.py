"""
Classifies: CHEBI:18179 phosphoinositide
"""
from rdkit import Chem

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    A phosphoinositide is a phosphatidylinositol phosphorylated at one or more hydroxy groups of inositol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a phosphoinositide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone part of phosphatidylinositol
    glycerol_pattern = Chem.MolFromSmarts("C(C(O)COC(=O))")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone part of phosphatidylinositol found"

    # Check for inositol ring with hydroxy groups
    inositol_ring_pattern = Chem.MolFromSmarts("C1([C@@H]([C@H]([C@@H]([C@@H]([C@H]1O)O)O)O)O)O")
    if not mol.HasSubstructMatch(inositol_ring_pattern):
        return False, "No inositol ring with proper hydroxy groups found"

    # Check for phosphate groups attached to inositol
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(O)=O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)

    # At least one phosphate group must be attached to the inositol ring
    inositol_phosphate_pattern = Chem.MolFromSmarts("[C@H](OP(O)(O)=O)[C@H]1([C@@H]([C@H]([C@@H]([C@@H]([C@H]1O)O)O)O)O)O")
    if not mol.HasSubstructMatch(inositol_phosphate_pattern) or len(phosphate_matches) < 1:
        return False, "No proper phosphate groups attached to inositol"

    return True, "Contains glycerol backbone, inositol ring, and phosphorylated groups specific to phosphoinositides"