"""
Classifies: CHEBI:18179 phosphoinositide
"""
"""
Classifies: CHEBI:28874 phosphoinositide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    A phosphoinositide is a phosphatidylinositol with one or more phosphate groups attached to the inositol ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphoinositide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for inositol ring (6-membered ring with multiple hydroxyl groups)
    inositol_pattern = Chem.MolFromSmarts("[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H]1O")
    inositol_matches = mol.GetSubstructMatches(inositol_pattern)
    if len(inositol_matches) == 0:
        return False, "No inositol ring found"

    # Look for phosphate groups attached to inositol
    inositol_phosphate_pattern = Chem.MolFromSmarts("[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H]1OP(=O)(O)O")
    inositol_phosphate_matches = mol.GetSubstructMatches(inositol_phosphate_pattern)
    if len(inositol_phosphate_matches) == 0:
        return False, "No phosphate groups attached to inositol"

    # Look for at least one ester group (indicative of a lipid tail)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) == 0:
        return False, "No ester groups found"

    return True, "Contains inositol ring with at least one phosphate group and ester groups"