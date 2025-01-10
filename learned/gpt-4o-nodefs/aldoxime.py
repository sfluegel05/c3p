"""
Classifies: CHEBI:22307 aldoxime
"""
from rdkit import Chem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime based on its SMILES string.
    An aldoxime contains a carbon-nitrogen double bond (C=N) and an oxygen 
    single-bonded to hydrogen (OH group) forming part of an aldehyde's oxime.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aldoxime, False otherwise
        str: Reason for classification
    """

    # Parse SMILES into RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for an aldoxime (R-CH=NOH)
    aldoxime_pattern = Chem.MolFromSmarts("[#6:1]=[N:2][O:3][H]")

    if mol.HasSubstructMatch(aldoxime_pattern):
        return True, "Contains the characteristic C=N-OH structure of aldoximes"
    
    return False, "Does not contain the characteristic C=N-OH structure of aldoximes"