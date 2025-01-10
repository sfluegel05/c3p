"""
Classifies: CHEBI:2440 acrovestone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_acrovestone(smiles: str):
    """
    Attempts to determine if a molecule is an acrovestone-like structure based on its SMILES string.
    We hypothesize that acrovestone involves isoflavone structures with sugar moieties.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule shows key features of acrovestone-like structures, False otherwise
        str: Reason for classification or lack thereof
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for isoflavone core structure (estimated 4H-chromen-4-one) presence
    isoflavone_pattern = Chem.MolFromSmarts("C1=CC(=CC=C1C=2C(=O)C3=CC=CC=C3O2)")
    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "No isoflavone core structure found"
    
    # Check for sugar moieties (glycosidic bonds)
    glycoside_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H]([C@@H](O)[C@H](O)[C@H]1O)")  # Example for glycoside
    if not mol.HasSubstructMatch(glycoside_pattern):
        return False, "No recognizable glycosidic bonds (sugar moieties) attached"
    
    return True, "Typical acrovestone-like structure detected: isoflavone with glycosidic linkage"