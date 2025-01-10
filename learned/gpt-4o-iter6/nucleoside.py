"""
Classifies: CHEBI:33838 nucleoside
"""
from rdkit import Chem

def is_nucleoside(smiles: str):
    """
    Determines if a molecule is a nucleoside based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a nucleoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Expand nucleobase patterns to include common tautomeric forms
    nucleobase_patterns = [
        Chem.MolFromSmarts("c1ncnc2ncnc12"),  # Purine base
        Chem.MolFromSmarts("c1[nH]c(=O)[nH]c(=O)[nH]1"),  # Uracil/Thymine
        Chem.MolFromSmarts("c1nc2c([nH]1)ncnc2=O"),  # Guanine
        Chem.MolFromSmarts("c1[nH]c(=O)ncn1"),  # Cytosine
        Chem.MolFromSmarts("c1nc2[nH]c[nH]c2[nH]1")  # Adenine
    ]
    
    # Check for nucleobase substructures
    if not any(mol.HasSubstructMatch(pattern) for pattern in nucleobase_patterns):
        return False, "No nucleobase found"
    
    # Ribose and deoxyribose sugar patterns
    ribose_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@@H]([C@H]([C@@H]1O)O)C")  # Ribose
    deoxyribose_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@H]([C@H]([C@H]1)O)C")  # Deoxyribose - correct stereo flexibility allowed

    # Check for sugar substructures
    if not (mol.HasSubstructMatch(ribose_pattern) or mol.HasSubstructMatch(deoxyribose_pattern)):
        return False, "No ribose or deoxyribose sugar found"
    
    # Successful identification if both a nucleobase and a ribose/deoxyribose are found
    return True, "Contains a nucleobase and a ribose or deoxyribose sugar with appropriate linkage"