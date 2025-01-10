"""
Classifies: CHEBI:18241 2'-deoxyribonucleoside 5'-monophosphate
"""
from rdkit import Chem

def is_2__deoxyribonucleoside_5__monophosphate(smiles: str):
    """
    Determines if a molecule is a 2'-deoxyribonucleoside 5'-monophosphate based on its SMILES string.
    These molecules have a 2'-deoxyribose sugar, a phosphate group at the 5' position,
    and a nucleobase attached to the 1' position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2'-deoxyribonucleoside 5'-monophosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (None, "Invalid SMILES string")
    
    # Check for 2'-deoxyribose moiety (correcting pattern and allowing flexibility)
    deoxyribose_pattern = Chem.MolFromSmarts('[C@H]1([C@@H](O)[C@H](CO)O)[C@@H]2O[C@H](C1)O2')
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return (False, "No 2'-deoxyribose moiety found")
    
    # Correct phosphate group pattern at 5' position
    phosphate_pattern = Chem.MolFromSmarts('O[P](=O)(O)O[C@H]1[C@H](O)[C@@H](O).C1')
    if not mol.HasSubstructMatch(phosphate_pattern):
        return (False, "No 5'-phosphate group found")
        
    # Broadened patterns for nucleobase detection (common structures)
    nucleobase_patterns = [
        Chem.MolFromSmarts('n1c2c([nH]c([nH]1)N)c([nH]c2)N'),  # adenine
        Chem.MolFromSmarts('C1=NC2=C(N[C@H]3O[C@H](COP(O)(=O)O)[C@@H](O)C3)N=C(N)N2CN1'),  # guanine
        Chem.MolFromSmarts('n1cc(Cc2cc[nH]c2)[nH]c(=O)c1'),    # cytosine
        Chem.MolFromSmarts('n1c(C=O)c[nH]c(=O)c1'),            # uracil and thymine
    ]

    nucleobase_found = any(mol.HasSubstructMatch(pattern) for pattern in nucleobase_patterns)
    
    if not nucleobase_found:
        return (False, "No suitable nucleobase found attached to the sugar")
    
    return (True, "Contains 2'-deoxyribose, 5'-phosphate group, and a nucleobase")