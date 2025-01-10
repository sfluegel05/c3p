"""
Classifies: CHEBI:2440 acrovestone
"""
from rdkit import Chem

def is_acrovestone(smiles: str):
    """
    Classifies if a molecule is related to acrovestone based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is considered an acrovestone-related structure, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Generic isoflavone-like ring core pattern
    isoflavone_core_pattern = Chem.MolFromSmarts("c1ccc2c(c1)C(=O)c3cc(O)ccc3O2")
    if not mol.HasSubstructMatch(isoflavone_core_pattern):
        return False, "No broad isoflavone-like core structure found"

    # Flexible glycosidic linkage pattern allowing for diverse sugars
    glycoside_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O1)C")
    if not mol.HasSubstructMatch(glycoside_pattern):
        return False, "No flexible glycosidic linkage found"
    
    # Check for the presence of hydroxy and/or methoxy groups flexibly
    substitutions_pattern_hydroxy = Chem.MolFromSmarts("[OX2H]")
    substitutions_pattern_methoxy = Chem.MolFromSmarts("C[OX2H]")
    
    has_hydroxy = mol.HasSubstructMatch(substitutions_pattern_hydroxy)
    has_methoxy = mol.HasSubstructMatch(substitutions_pattern_methoxy)
    
    if not (has_hydroxy or has_methoxy):
        return False, "No relevant hydroxy or methoxy substitutions found"
    
    return True, "Matches acrovestone structure with flexible isoflavone core and polyphenolic substituents"