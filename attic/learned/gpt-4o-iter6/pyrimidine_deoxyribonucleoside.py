"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
from rdkit import Chem

def is_pyrimidine_deoxyribonucleoside(smiles: str):
    """
    Determines if a molecule is a pyrimidine deoxyribonucleoside based on its SMILES string.
    A pyrimidine deoxyribonucleoside is defined as a deoxynucleoside containing a pyrimidine base.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrimidine deoxyribonucleoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Deoxyribose (2'-deoxyribose) pattern - common sugar in DNA
    deoxyribose_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@H](CO)[C@H]1O")
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No deoxyribose sugar found"
    
    # Pyrimidine base pattern - includes uracil, thymine, or cytosine
    pyrimidine_pattern_uracil = Chem.MolFromSmarts("C1=NC(=O)NC(=O)C=C1")
    pyrimidine_pattern_thymine = Chem.MolFromSmarts("C1=NC(=O)N(C)C(=O)C=C1")
    pyrimidine_pattern_cytosine = Chem.MolFromSmarts("C1=NC(=O)NC(N)=C1")
    
    if not (mol.HasSubstructMatch(pyrimidine_pattern_uracil) or
            mol.HasSubstructMatch(pyrimidine_pattern_thymine) or
            mol.HasSubstructMatch(pyrimidine_pattern_cytosine)):
        return False, "No pyrimidine base found"

    return True, "Contains both deoxyribose sugar and pyrimidine base"