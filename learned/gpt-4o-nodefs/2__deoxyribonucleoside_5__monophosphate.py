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
    
    # Check for 2'-deoxyribose moiety (allowing variable stereochemistry)
    deoxyribose_pattern = Chem.MolFromSmarts('[C@H]1C[C@@H](O)[C@H](CO)O1')
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return (False, "No 2'-deoxyribose moiety found")

    # Check for a phosphate group at 5' position
    phosphate_pattern = Chem.MolFromSmarts('COP(=O)(O)O')
    if not mol.HasSubstructMatch(phosphate_pattern):
        return (False, "No 5'-phosphate group found")
        
    # Check for nucleobase attached to the sugar via a variable pattern
    nucleobase_pattern = Chem.MolFromSmarts('n')
    if not mol.HasSubstructMatch(nucleobase_pattern):
        return (False, "No nucleobase found (purine/pyrimidine) attached to the sugar")
    
    return (True, "Contains 2'-deoxyribose, 5'-phosphate group, and a nucleobase")