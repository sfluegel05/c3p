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
    
    # Check for 2'-deoxyribose moiety 
    # A more flexible pattern allowing both alpha and beta anomers and specific deoxyribose configuration
    deoxyribose_pattern = Chem.MolFromSmarts('[C@@H]1(O)[C@H](O)[C@@H](CO)O[C@H](C1)O')
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return (False, "No 2'-deoxyribose moiety found")
    
    # Check for a phosphate group at 5' position (more anchoring to deoxyribose hydroxyl)
    phosphate_pattern = Chem.MolFromSmarts('[O-]P(=O)(O)O[C@H])')  # indicates 5' terminal position 
    if not mol.HasSubstructMatch(phosphate_pattern):
        return (False, "No 5'-phosphate group found")
        
    # Check for nucleobase attached to the sugar
    # More specific nucleobase patterns - adenine, guanine, cytosine, thymine, and uracil
    nucleobase_patterns = [Chem.MolFromSmarts('n1cnc2c1ncnc2'),  # adenine
                           Chem.MolFromSmarts('n1cnc2c1[nH]cn2'),  # guanine
                           Chem.MolFromSmarts('n1ccn(C)c(=O)c1'),   # cytosine
                           Chem.MolFromSmarts('c1ncnc2[nH]ccc12'),  # uracil
                           Chem.MolFromSmarts('c1ccn(C)c(=O)n1')]   # thymine

    nucleobase_found = any(mol.HasSubstructMatch(pattern) for pattern in nucleobase_patterns)
    
    if not nucleobase_found:
        return (False, "No nucleobase found (purine/pyrimidine) attached to the sugar")
    
    return (True, "Contains 2'-deoxyribose, 5'-phosphate group, and a nucleobase")