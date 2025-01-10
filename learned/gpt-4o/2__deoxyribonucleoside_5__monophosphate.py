"""
Classifies: CHEBI:18241 2'-deoxyribonucleoside 5'-monophosphate
"""
from rdkit import Chem

def is_2__deoxyribonucleoside_5__monophosphate(smiles: str):
    """
    Determines if a molecule is a 2'-deoxyribonucleoside 5'-monophosphate based on its SMILES string.
    A 2'-deoxyribonucleoside monophosphate contains a deoxyribose sugar, a nitrogenous base, and a
    phosphate group at the 5' position.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule fits the class, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for deoxyribose sugar pattern (remove hydroxyl group at 2' position)
    # Adjusted pattern to be more flexible in atom valence and stereochemistry
    deoxyribose_pattern = Chem.MolFromSmarts("C1OC(CO)C(O)C1")
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No 2'-deoxyribose structure found"

    # Check for phosphate group pattern at 5' position
    # Broadened pattern to match potential variations in bonding
    phosphate_pattern = Chem.MolFromSmarts("COP(=O)(O)[O-]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        phosphate_pattern = Chem.MolFromSmarts("COP(=O)(O)O")
        if not mol.HasSubstructMatch(phosphate_pattern):
            return False, "No 5'-phosphate group found"

    # Check for typical nucleobases pattern - more permutations
    base_patterns = [
        Chem.MolFromSmarts("n1cnc2c1ncnc2"),  # adenine base
        Chem.MolFromSmarts("n1c(O)nc2c1ncnc2"),  # guanine base
        Chem.MolFromSmarts("n1cc(nc1)C=O"),  # cytosine base
        Chem.MolFromSmarts("c1ccn(c1)C=O"),  # thymine base
        Chem.MolFromSmarts("n1[CH]cc(=O)[nh]c1=O")  # uracil base
    ]
    if not any(mol.HasSubstructMatch(base) for base in base_patterns):
        return False, "No typical nucleobase structure found"

    return True, "Structure fits 2'-deoxyribonucleoside 5'-monophosphate definition"