"""
Classifies: CHEBI:36976 nucleotide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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
    # Comprehensive pattern for phosphate group (consider multiple esterifications)
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)"),
    esterified_phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")

    # Patterns for ribose or deoxyribose sugar backbone
    sugar_ribose_pattern = Chem.MolFromSmarts("C1(O)C(O)C(O)C(O1)")
    sugar_deoxyribose_pattern = Chem.MolFromSmarts("C1(O)C(O)C(CO)O1")

    # Patterns for nitrogenous bases (including variations in purines and pyrimidines)
    base_purine_pattern = Chem.MolFromSmarts("n1cnc2ncnc2c1")
    base_pyrimidine_pattern = Chem.MolFromSmarts("n1cc(C=O)cn1")

    # Check for presence of phosphate group
    if not (mol.HasSubstructMatch(phosphate_pattern) or mol.HasSubstructMatch(esterified_phosphate_pattern)):
        return False, "No phosphate group found"

    # Check for presence of sugar backbone (consider both ribose and deoxyribose)
    if not (mol.HasSubstructMatch(sugar_ribose_pattern) or mol.HasSubstructMatch(sugar_deoxyribose_pattern)):
        return False, "No sugar backbone found"
    
    # Check for presence of a nitrogenous base
    if not (mol.HasSubstructMatch(base_purine_pattern) or mol.HasSubstructMatch(base_pyrimidine_pattern)):
        return False, "No recognizable nitrogenous base found"

    return True, "Contains a sugar backbone, nitrogenous base, and phosphate group"