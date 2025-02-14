"""
Classifies: CHEBI:46633 carbapenems
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_carbapenems(smiles: str):
    """
    Determines if a molecule is a carbapenem based on its SMILES string.
    Carbapenems are a class of beta-lactam antibiotics with a specific bicyclic structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbapenem, False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core carbapenem structure using SMARTS
    # Allow for variation in the N position and ring double bond and sterochemistry
    carbapenem_core_smarts1 = "[N]1[C@H]2[C@H](C1=O)C[C@@H]2"
    carbapenem_core_smarts2 = "[C]1[C@H]2[C@H](N1C=O)C[C@@H]2"


    carbapenem_core_pattern1 = Chem.MolFromSmarts(carbapenem_core_smarts1)
    carbapenem_core_pattern2 = Chem.MolFromSmarts(carbapenem_core_smarts2)


    if not (mol.HasSubstructMatch(carbapenem_core_pattern1) or mol.HasSubstructMatch(carbapenem_core_pattern2)):
       return False, "Core carbapenem structure not found"

    # Check for a carboxylic acid group (C(=O)O), which most carbapenems have
    acid_group_smarts = "C(=O)O"
    acid_group_pattern = Chem.MolFromSmarts(acid_group_smarts)

    if not mol.HasSubstructMatch(acid_group_pattern):
       return False, "No carboxylic acid group present"


    return True, "Contains the core carbapenem structure and an acid group"