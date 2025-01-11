"""
Classifies: CHEBI:36976 nucleotide
"""
from rdkit import Chem

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.
    A nucleotide generally has a sugar backbone, a nitrogenous base, and one or more phosphate groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleotide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define patterns for nucleotide detection
    # Comprehensive pattern for phosphate group (potential for multiple links and esterification)
    phosphate_pattern = Chem.MolFromSmarts("P(O)(=O)(O)")

    # Patterns for ribose or deoxyribose sugar backbone (considering stereochemistry)
    sugar_ribose_pattern = Chem.MolFromSmarts("[C@H]1([C@H](O)[C@@H](O)[C@H](O1)CO)")
    sugar_deoxyribose_pattern = Chem.MolFromSmarts("[C@H]1([C@H](O)[C@@H](O)[C@H](C1)CO)")

    # Patterns for nitrogenous bases (including variations)
    purine_base_pattern = Chem.MolFromSmarts("n1cnc2ncnc2c1")
    pyrimidine_base_pattern = Chem.MolFromSmarts("n1ccnc1")

    # Check for presence of phosphate group
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group found"
    
    # Check for presence of sugar backbone (ribose or deoxyribose)
    sugar_matches = mol.GetSubstructMatches(sugar_ribose_pattern) or mol.GetSubstructMatches(sugar_deoxyribose_pattern)
    if not sugar_matches:
        return False, "No sugar backbone found"
    
    # Check for presence of a nitrogenous base (purines or pyrimidines)
    base_matches = mol.GetSubstructMatches(purine_base_pattern) or mol.GetSubstructMatches(pyrimidine_base_pattern)
    if not base_matches:
        return False, "No recognizable nitrogenous base found"

    return True, "Contains a sugar backbone, nitrogenous base, and phosphate group"