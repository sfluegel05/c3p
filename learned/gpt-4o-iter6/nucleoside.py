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
    
    # Nucleobase patterns (expanded coverage)
    purine_patterns = [
        Chem.MolFromSmarts("c1[nH]cnc2ncnc12"), # Generic purine
        Chem.MolFromSmarts("n1cnc2c([nH]1)ncnc2"), # Adenine-like
        Chem.MolFromSmarts("n1cnc2c(nc[nH]2)[nH]1") # Guanine-like
    ]
    
    pyrimidine_patterns = [
        Chem.MolFromSmarts("c1c[nH]c(=O)[nH]c1"), # Generic pyrimidine
        Chem.MolFromSmarts("c1cnc(=O)[nH]c1"), # Uracil-like
        Chem.MolFromSmarts("c1c[nH]cnc1=O") # Thymine-like
    ]
    
    # Check for nucleobase substructures
    has_nucleobase = any(mol.HasSubstructMatch(pattern) for pattern in purine_patterns + pyrimidine_patterns)
    
    if not has_nucleobase:
        return False, "No nucleobase found"
    
    # Ribose and deoxyribose sugar patterns including stereo
    ribose_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@H](CO)O[C@@H]1") # Ribose with stereo
    deoxyribose_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@H](C)O[C@@H]1") # Deoxyribose with stereo
    
    # Check for sugar substructures
    has_sugar = mol.HasSubstructMatch(ribose_pattern) or mol.HasSubstructMatch(deoxyribose_pattern)
    
    if not has_sugar:
        return False, "No ribose or deoxyribose sugar found"
    
    # Check N-glycosidic bond
    # Generalized pattern for glycosidic bond with nitrogen atom
    n_glycosidic_bond = Chem.MolFromSmarts("[O][C@H]1[C@@H]([C@H]([C@H](O1)CO)O)N")
    has_glycosidic_bond = mol.HasSubstructMatch(n_glycosidic_bond)
    
    if not has_glycosidic_bond:
        return False, "Missing N-glycosidic linkage between sugar and nucleobase"
    
    return True, "Contains a nucleobase and a ribose or deoxyribose sugar with appropriate linkage"

# Example SMILES testing, should test against a list of known nucleoside SMILES for validation