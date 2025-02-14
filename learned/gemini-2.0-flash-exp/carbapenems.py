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
    # This pattern correctly identifies the fused ring system and the beta-lactam
    carbapenem_core_smarts = "[C]1[C]([C](=[O])[N]2[C]1[C][C]2)"
    carbapenem_core_pattern = Chem.MolFromSmarts(carbapenem_core_smarts)

    if not mol.HasSubstructMatch(carbapenem_core_pattern):
       return False, "Core carbapenem structure not found"
    
    return True, "Contains the core carbapenem structure"