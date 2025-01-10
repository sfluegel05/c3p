"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate based on its SMILES string.
    Phosphatidylinositol phosphates are characterized by a glycerol backbone with two fatty acid chains,
    a phosphate group, and an inositol ring that is phosphorylated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is phosphatidylinositol phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern
    glycerol_pattern = Chem.MolFromSmarts("[O]C[C@H](CO[H])")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
        
    # Look for fatty acid ester groups
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C@H]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"

    # Look for inositol ring with phosphate
    inositol_phosphate_pattern = Chem.MolFromSmarts("C1C(O)C(O)C(O)C(O)C1OP(=O)(O)O")
    if not mol.HasSubstructMatch(inositol_phosphate_pattern):
        return False, "Missing inositol ring with phosphate group"

    # Check for phosphate groups (should have at least one)
    phosphate_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("OP(=O)(O)O")))
    if phosphate_count < 1:
        return False, "Insufficient number of phosphate groups"

    return True, "Contains features characteristic of phosphatidylinositol phosphate"