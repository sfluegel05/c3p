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
    
    # Isoflavone core structure with room for substitution allowances
    isoflavone_pattern = Chem.MolFromSmarts("c1cc2oc3c(ccc4ccccc34)c(=O)c2c[cH]1")
    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "No isoflavone core structure found"

    # Simplified glycosidic linkage pattern/common sugar units
    glycoside_pattern = Chem.MolFromSmarts("[OX2H,[OX2CH3]R]1[C@H]([C@@H](O)C(O)C1O[*])")
    if not mol.HasSubstructMatch(glycoside_pattern):
        return False, "No glycosidic linkage found"
    
    # Flexible hydroxy and methoxy group patterns on aromatic rings
    hydroxy_methoxy_pattern = Chem.MolFromSmarts("[c][OX2H] | [c][OX2CH3]")
    if not mol.HasSubstructMatch(hydroxy_methoxy_pattern):
        return False, "No relevant hydroxy or methoxy substitutions found"
    
    return True, "Matches isoflavone core with glycosidic linkage and typical polyphenol substituents"