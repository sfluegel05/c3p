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
    
    # Look for core isoflavone structure
    isoflavone_pattern = Chem.MolFromSmarts("c1cc2c(cc1)C(=O)c3ccccc3O2")
    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "No isoflavone core structure found"
    
    # Check for glycosidic linkages
    glycoside_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](O[C@@H]1[*])CO")
    if not mol.HasSubstructMatch(glycoside_pattern):
        return False, "No glycosidic linkage found"
    
    # Look for hydroxy and methoxy groups on aromatic rings
    hydroxy_methoxy_pattern = Chem.MolFromSmarts("[OX2H][cR1] | [OX2CH3][cR1]")
    if not mol.HasSubstructMatch(hydroxy_methoxy_pattern):
        return False, "No relevant hydroxy or methoxy substitutions found"
    
    return True, "Matches isoflavone core with glycosidic linkage and typical polyphenol substituents"