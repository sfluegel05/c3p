"""
Classifies: CHEBI:72600 spiroketal
"""
from rdkit import Chem

def is_spiroketal(smiles: str):
    """
    Determines if the given SMILES string represents a spiroketal.
    
    A spiroketal is characterized by:
    - A spiro junction connecting at least two rings.
    - A cyclic ketal group as part of the compound.
    
    Args:
        smiles (str): The SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a spiroketal, False otherwise.
        str: Explanation of the classification.
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"
    
    # Define complex spiroketal detection logic using SMARTS
    # This pattern accounts for spiro connection and ketal presence in a single molecule scope
    spiroketal_pattern = '[$([C]@/1([C])([O])[C][C][O]1),$([C]@[C]1[C][O][C][C]1)]'
    
    if not mol.HasSubstructMatch(Chem.MolFromSmarts(spiroketal_pattern)):
        return False, "No spiroketal features found in the molecule"
    
    return True, "Contains a spiro junction with a cyclic ketal group"