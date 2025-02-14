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
    # Allow for variation in the N position, bond order, and any substitution
    carbapenem_core_smarts = "[N,C]1[C]2[C](C1=O)[C]([C]=C2)"
    carbapenem_core_pattern = Chem.MolFromSmarts(carbapenem_core_smarts)

    if not mol.HasSubstructMatch(carbapenem_core_pattern):
       return False, "Core carbapenem structure not found"

    # Check for a carboxylic acid or carboxylate group (C(=O)[O,O-]), which most carbapenems have
    # Make this optional as some can be pro-drugs or sodium salts
    acid_group_smarts = "C(=O)[O,O-]"
    acid_group_pattern = Chem.MolFromSmarts(acid_group_smarts)

    #if not mol.HasSubstructMatch(acid_group_pattern):
    #   return False, "No carboxylic acid or carboxylate group present"
    
    return True, "Contains the core carbapenem structure"