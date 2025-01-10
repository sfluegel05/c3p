"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
from rdkit import Chem

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
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")  # Basic glycerol pattern
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
        
    # Recognize two fatty acid ester linkages on glycerol
    ester_pattern = Chem.MolFromSmarts("C(=O)OCC")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester linkages, need at least 2"

    # Inositol ring identification with variable positions of phosphates
    inositol_pattern = Chem.MolFromSmarts("C1(O)C(O)C(O)C(O)C(O)C1OP")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "Inositol ring with phosphates not found"

    # Check phosphate groups on inositol ring
    phosphate_groups_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    phosphate_count = len(mol.GetSubstructMatches(phosphate_groups_pattern))
    if phosphate_count < 1:
        return False, "Insufficient number of phosphate groups on inositol"

    return True, "Structure consistent with phosphatidylinositol phosphate"