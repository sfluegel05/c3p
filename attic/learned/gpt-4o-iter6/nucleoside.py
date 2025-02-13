"""
Classifies: CHEBI:33838 nucleoside
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

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
    
    # Nucleobase patterns (adenine, guanine, cytosine, uracil, thymine)
    nucleobase_pattern = Chem.MolFromSmarts("n1cnc2c1ncnc2") # purine base pattern (e.g., adenine, guanine)
    pyrimidine_pattern = Chem.MolFromSmarts("c1c[nH]c(=O)n(c1=O)") # pyrimidine base pattern (e.g., cytosine, uracil, thymine)
    
    has_purine_like = mol.HasSubstructMatch(nucleobase_pattern)
    has_pyrimidine_like = mol.HasSubstructMatch(pyrimidine_pattern)
    
    if not has_purine_like and not has_pyrimidine_like:
        return False, "No nucleobase found"
    
    # Sugar pattern (ribose or deoxyribose)
    sugar_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@H](O)[C@@H](CO)C1")
    deoxyribose_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@H](O)[C@@H](C)C1") # for deoxyribose
    
    has_sugar = mol.HasSubstructMatch(sugar_pattern) or mol.HasSubstructMatch(deoxyribose_pattern)
    
    if not has_sugar:
        return False, "No ribose or deoxyribose sugar found"
    
    return True, "Contains both a nucleobase and a ribose or deoxyribose sugar"

# Example SMILES testing:
# print(is_nucleoside("O=C1N(C=CC(N1)=O)C2OC(CO)C(C2O)O"))  # Expected True, Spongouridin