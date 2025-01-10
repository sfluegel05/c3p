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
    
    # Broadened polyphenolic core pattern to include diverse isoflavone derivatives
    polyphenolic_core_pattern = Chem.MolFromSmarts("c1cc(O)c2c(c1)c(=O)c3c(O)cc(O)cc3o2")
    if not mol.HasSubstructMatch(polyphenolic_core_pattern):
        return False, "No representative polyphenolic core structure found"

    # Broadened pattern for glycoside linkage with diverse sugars
    glycoside_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O1)C")
    if not mol.HasSubstructMatch(glycoside_pattern):
        return False, "No flexible glycosidic linkage found"
    
    # Check for the presence of hydroxy and/or methoxy groups on an aromatic ring
    hydroxy_methoxy_pattern = Chem.MolFromSmarts("[OX2H]c|C[OX2H]c")
    
    if not mol.HasSubstructMatch(hydroxy_methoxy_pattern):
        return False, "No relevant hydroxy or methoxy aromatic substitutions found"
    
    return True, "Matches acrovestone structure with broadened polyphenolic core and diverse glycosidic linkages"