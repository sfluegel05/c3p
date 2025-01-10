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
    
    # Generalized isoflavone core structure pattern
    # Allowing for flexibility for substitution
    isoflavone_pattern = Chem.MolFromSmarts("c1ccc2[o,c]c(=O)[c,C]c3c(c2c1)c(=O)c[c,C]4[c,C][c,C][c,C][c,C][c,C]4")
    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "No isoflavone core structure found"

    # Flexible glycosidic linkage or sugar moiety patterns
    # Accommodating common sugar units and linkages
    glycoside_pattern = Chem.MolFromSmarts("O[C@@H]1[C@H]([C@H]([C@@H](O)[C@@H]1O)O[*])")
    if not mol.HasSubstructMatch(glycoside_pattern):
        return False, "No suitable glycosidic linkage found"
    
    # Flexible check for hydroxy and methoxy groups
    hydroxyl_methoxy_pattern = Chem.MolFromSmarts("[c][OH] | [c][OMe] | [c][OX2H] | [c][OX2CH3]")
    if not mol.HasSubstructMatch(hydroxyl_methoxy_pattern):
        return False, "No relevant hydroxy or methoxy substitutions found"

    return True, "Matches isoflavone core with glycosidic linkage and typical polyphenol substituents"