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
    # Pattern for a phosphate group (pat => P(=O)(O)O)
    phosphate_pattern = Chem.MolFromSmarts("[$(P(=O)([O-])([O-])[O-]), $(P(=O)(O)(O))]")

    # Pattern for a ribose or deoxyribose sugar backbone (H bonded to C in a chain)
    sugar_pattern = Chem.MolFromSmarts("C1(O)([CH2])OC[C@H]1O")
    
    # Pattern for a nitrogenous base (detecting purines/pyrimidines is complex)
    base_pattern = Chem.MolFromSmarts("n1cnc2c1ncnc2")

    # Check for presence of phosphate group
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group found"
    
    # Check for presence of sugar backbone
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No sugar backbone found"
    
    # Check for presence of nitrogenous base
    base_matches = mol.GetSubstructMatches(base_pattern)
    if not base_matches:
        return False, "No recognizable nitrogenous base found"

    return True, "Contains a sugar backbone, nitrogenous base, and phosphate group"