"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
from rdkit import Chem

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.
    A nucleoside 5'-phosphate consists of a ribose or deoxyribose sugar linked to a purine or pyrimidine base,
    with phosphorylation at the C-5 position of the sugar. 

    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a nucleoside 5'-phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define improved sugar patterns for ribose and deoxyribose including stereochemistry
    ribose_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@@H]([C@@H](O)C1)O")
    deoxyribose_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@@H]([C@@H](C)O1)")

    # Improved phosphate group pattern for multiple phosphorylations at C-5 position
    phosphate_pattern = Chem.MolFromSmarts("[O][P](=O)([O])[O][CH2]")

    # Revised patterns for purine and pyrimidine bases allowing for common derivatives
    purine_pattern = Chem.MolFromSmarts("n1cnc2c1ncnc2")
    pyrimidine_pattern = Chem.MolFromSmarts("n1cncc1=O")

    # Check for sugar pattern matches in the molecule
    if not mol.HasSubstructMatch(ribose_pattern) and not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "Does not contain a ribose or deoxyribose sugar structure"

    # Check for phosphate group patterns
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found attached at 5' position of the sugar"

    # Check for either purine or pyrimidine base
    if not mol.HasSubstructMatch(purine_pattern) and not mol.HasSubstructMatch(pyrimidine_pattern):
        return False, "No recognized purine or pyrimidine base structure found"

    return True, "Contains a ribosyl or deoxyribosyl sugar, phosphate at 5', and a nucleobase"