"""
Classifies: CHEBI:47622 acetate ester
"""
"""
Classifies: CHEBI:48030 acetate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_acetate_ester(smiles: str):
    """
    Determines if a molecule is an acetate ester based on its SMILES string.
    An acetate ester is any carboxylic ester where the carboxylic acid component is acetic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acetate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the acetate ester pattern: -O-C(=O)-CH3
    acetate_ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[CX4H3]")
    
    # Check if the pattern is present in the molecule
    if mol.HasSubstructMatch(acetate_ester_pattern):
        return True, "Contains the acetate ester functional group (-O-C(=O)-CH3)"
    else:
        return False, "No acetate ester functional group found"

# Example usage:
# print(is_acetate_ester("CCCCOC(C)=O"))  # Should return (True, "Contains the acetate ester functional group (-O-C(=O)-CH3)")