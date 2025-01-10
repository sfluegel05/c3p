"""
Classifies: CHEBI:35915 sterol ester
"""
from rdkit import Chem

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    A steroid ester is obtained by formal condensation of 
    the carboxy group of any carboxylic acid with the 3-hydroxy group of a sterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a general sterol backbone pattern (cyclopentanoperhydrophenanthrene)
    sterol_pattern = Chem.MolFromSmarts("C1CCC2=C(C1)C3CCC4C(C=CC(O)C4)(C3C2)C")
    if not mol.HasSubstructMatch(sterol_pattern):
        return False, "No sterol backbone found"

    # Define ester linkage pattern (ester with -OH on sterol)
    ester_linkage_pattern = Chem.MolFromSmarts("C(=O)O[C][C][C]O")
    if not mol.HasSubstructMatch(ester_linkage_pattern):
        return False, "No appropriate ester linkage found"

    return True, "Contains sterol backbone with appropriate ester linkage"