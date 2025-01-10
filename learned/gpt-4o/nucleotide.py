"""
Classifies: CHEBI:36976 nucleotide
"""
from rdkit import Chem

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.
    A nucleotide is a nucleoside phosphate resulting from the condensation of the 3 or 5 hydroxy group of a nucleoside with phosphoric acid.

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

    # Define SMARTS patterns for nucleobases (purines and pyrimidines), sugar, and phosphate groups
    purine_pattern = Chem.MolFromSmarts("n1cnc2ncnc12")  # Recognize purine rings
    pyrimidine_pattern = Chem.MolFromSmarts("n1ccc(N)nc1=O")  # Recognize pyrimidine rings
    # Updated sugar pattern to match ribose and deoxyribose more generally
    sugar_pattern = Chem.MolFromSmarts("[C@H]1([C@H]([C@H]([C@H](O1)CO)O)O)O")  # All parts of ribose/deoxyribose with stereo variations
    # Updated phosphate pattern to consider multiple phosphates with various bond/orders
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)([O-])[O-]")  # Recognize basic phosphate structure

    # Check for nucleobase presence (only basic check of purine or pyrimidine)
    if not (mol.HasSubstructMatch(purine_pattern) or mol.HasSubstructMatch(pyrimidine_pattern)):
        return False, "No nucleobase found"

    # Check for sugar structure presence
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar structure (ribose or deoxyribose) found"
    
    # Check for phosphate group presence
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    return True, "Contains nucleobase attached to a sugar with a phosphate group"