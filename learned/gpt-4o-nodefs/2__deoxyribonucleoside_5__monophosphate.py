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
        return False, "Invalid SMILES string"
    
    # Check for 2'-deoxyribose moiety (ensuring correct stereochemistry)
    deoxyribose_pattern = Chem.MolFromSmarts('[C@@H]1[C@H]([C@H](C(O1)O)O)C')
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No 2'-deoxyribose moiety found"

    # Check for a phosphate group at 5' position
    phosphate_pattern = Chem.MolFromSmarts('COP(=O)(O)O')
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No 5'-phosphate group found"
        
    # Check for nucleobase attached to the anomeric carbon (C1')
    # This can be various nucleobases, represented by generic aromatic nitrogen-containing heterocycles
    nucleobase_pattern = Chem.MolFromSmarts('n1cnc[nH]c1')  # Generic purine/pyrimidine ring structure
    if not mol.HasSubstructMatch(nucleobase_pattern):
        return False, "No nucleobase (purine/pyrimidine) found attached to the anomeric carbon"

    # All checks passed
    return True, "Contains 2'-deoxyribose, 5'-phosphate group, and a nucleobase"