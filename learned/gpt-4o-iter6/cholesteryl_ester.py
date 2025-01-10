"""
Classifies: CHEBI:17002 cholesteryl ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_cholesteryl_ester(smiles: str):
    """
    Determines if a molecule is a cholesteryl ester based on its SMILES string.
    A cholesteryl ester is defined by the ester linkage of any carboxylic acid with
    the 3-hydroxy group of cholesterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cholesteryl ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")
    
    # Improved cholesterol pattern
    cholesterol_pattern = Chem.MolFromSmarts("C1CCC2C1CCC3C2CCC4C3(CCCC4)C")
    if not mol.HasSubstructMatch(cholesterol_pattern):
        return (False, "No cholesterol framework identified")
    
    # Check for an ester group with connection to cholesterol
    ester_position_pattern = Chem.MolFromSmarts("O[C@H]1CC[C@H]2C[C@H]3CC[C@H]4[C@](CCC4)(CC3)C2CC1")
    if not mol.HasSubstructMatch(ester_position_pattern):
        return (False, "Ester linkage not correctly connected to the 3β-position of cholesterol")
    
    return (True, "SMILES indicates a cholesteryl ester structure")