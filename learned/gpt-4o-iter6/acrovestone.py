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
    
    # Isoflavone core with phenolic ring exemplified with C-C-C-C-C-O pattern
    # Used a more generic pattern to allow variability and substitution
    isoflavone_pattern = Chem.MolFromSmarts("C1=CC=C2C(=C1)C(=O)CC3=C2OCC3")
    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "No isoflavone core structure found"

    # Look for glycosidic linkage pattern (allow for diverse sugars, assume a simple ring with O-atoms)
    glycoside_pattern = Chem.MolFromSmarts("[OX2H][C@H]1[C@H]([C@H](O)[C@@H](O)[C@H]1O)")
    if not mol.HasSubstructMatch(glycoside_pattern):
        return False, "No common glycosidic linkage found"
    
    # Check for the presence of hydroxy groups
    hydroxy_pattern = Chem.MolFromSmarts("[CX3]=[CX3][OH]")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No relevant hydroxy substitutions found"

    # Check for the presence of methoxy groups
    methoxy_pattern = Chem.MolFromSmarts("[CX3]=[CX3][OCH3]")
    if not mol.HasSubstructMatch(methoxy_pattern):
        return False, "No relevant methoxy substitutions found"

    return True, "Matches isoflavone core with glycosidic linkage and common polyphenolic substituents"