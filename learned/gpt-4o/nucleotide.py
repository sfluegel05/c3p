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
        return None, "Invalid SMILES string"

    # Define SMARTS patterns for nucleobases (purines and pyrimidines)
    purine_pattern = Chem.MolFromSmarts("c1ncnc2[nH]c[nH]c12")  # General purine ring
    pyrimidine_pattern = Chem.MolFromSmarts("c1ncnc(N)c1=O")  # General pyrimidine ring

    # Define pattern for ribose or deoxyribose sugars, allowing for some variation
    sugar_pattern = Chem.MolFromSmarts("[C@@H]1(O)[C@@H](O)[C@H](CO)O[C@H]1")  # Ribose structure generalized

    # General phosphate pattern, accounting for potential variations
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)")

    # Check for the presence of nucleobase (purine or pyrimidine)
    if not (mol.HasSubstructMatch(purine_pattern) or mol.HasSubstructMatch(pyrimidine_pattern)):
        return False, "No nucleobase found"

    # Check for the presence of a sugar structure
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar structure (ribose or deoxyribose) found"
    
    # Check for the presence of at least one phosphate group
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    return True, "Molecule contains nucleobase, sugar, and phosphate group consistent with a nucleotide"