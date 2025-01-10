"""
Classifies: CHEBI:18241 2'-deoxyribonucleoside 5'-monophosphate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2__deoxyribonucleoside_5__monophosphate(smiles: str):
    """
    Determines if a molecule is a 2'-deoxyribonucleoside 5'-monophosphate based on its SMILES string.
    This class involves a 2'-deoxyribose sugar with a 5'-phosphate and one of the typical nucleobases.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2'-deoxyribonucleoside 5'-monophosphate, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Relaxed 2'-deoxyribose sugar pattern, focused on structural ring & connections
    deoxyribose_pattern = Chem.MolFromSmarts("C1OC(CO)CO1")
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No 2'-deoxyribose sugar structure found"

    # Phosphate group pattern with flexibility on charge
    phosphate_pattern = Chem.MolFromSmarts("COP(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No 5'-monophosphate group found"

    # Checking broadly for nucleobase structures
    # Supporting purine and pyrimidine base presence with simplified patterns
    purine_pattern = Chem.MolFromSmarts("n1cnc2c1ncnc2")
    pyrimidine_pattern = Chem.MolFromSmarts("c1ccn(c(=O)n1)")
    if not (mol.HasSubstructMatch(purine_pattern) or mol.HasSubstructMatch(pyrimidine_pattern)):
        return False, "No recognized purine or pyrimidine nucleobase found"

    return True, "Valid 2'-deoxyribonucleoside 5'-monophosphate with correct sugar, phosphate, and nucleobase"