"""
Classifies: CHEBI:18179 phosphoinositide
"""
from rdkit import Chem

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    A phosphoinositide is a phosphatidylinositol that is phosphorylated at one or more of the hydroxy groups of inositol.

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

    # Look for inositol ring pattern (6-membered ring with hydroxy groups)
    inositol_pattern = Chem.MolFromSmarts("C1C(O)C(O)C(O)C(O)C(O)C1")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"

    # Check for phosphate groups connected to inositol
    phosphate_pattern = Chem.MolFromSmarts("O[P](=O)([O-])[O-]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 1:
        return False, "No phosphorylation found on inositol"

    # Count the C(=O)O side chains indicating lipid linkage
    lipid_linkage_pattern = Chem.MolFromSmarts("C(=O)O[C@@H]")
    acyl_matches = mol.GetSubstructMatches(lipid_linkage_pattern)
    if not acyl_matches:
        return False, "No lipid linkage found"

    return True, "Contains inositol ring with phosphorylations and lipid linkages indicating a phosphoinositide"