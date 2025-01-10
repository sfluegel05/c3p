"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
from rdkit import Chem

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.
    Identifies purine or pyrimidine base, ribosyl or deoxyribosyl sugar, and phosphate group at the 5' position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside 5'-phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # General pattern for ribose and 2'-deoxyribose with more flexibility
    ribose_pattern = Chem.MolFromSmarts("OC1COC(CO)C1O")  # Looser match for ribose
    deoxyribose_pattern = Chem.MolFromSmarts("OC1COC(C)C1O") # Looser match for deoxyribose

    has_sugar = mol.HasSubstructMatch(ribose_pattern) or mol.HasSubstructMatch(deoxyribose_pattern)
    
    if not has_sugar:
        return False, "No ribose or deoxyribose sugar found"

    # More general pattern for Purine and Pyrimidine base to cover variations
    purine_pattern = Chem.MolFromSmarts("n1cnc2nc[nH]c2n1")  # A typical purine pattern
    pyrimidine_pattern = Chem.MolFromSmarts("n1c[nH]cnc1")  # A typical pyrimidine pattern

    has_base = mol.HasSubstructMatch(purine_pattern) or mol.HasSubstructMatch(pyrimidine_pattern)
    
    if not has_base:
        return False, "No purine or pyrimidine base found"

    # Phosphate group patterns with flexibility in position
    phosphate_pattern = Chem.MolFromSmarts("COP(O)(=O)O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    
    if not phosphate_matches:
        return False, "No phosphate group found at C-5'"

    return True, "Contains ribosyl or deoxyribosyl derivative of a purine/pyrimidine base with a 5'-phosphate group"