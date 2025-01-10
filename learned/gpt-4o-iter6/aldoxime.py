"""
Classifies: CHEBI:22307 aldoxime
"""
from rdkit import Chem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime based on its SMILES string.
    An aldoxime is defined as an oxime derived from an aldehyde, characterized by the structure RCH=NOH.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldoxime, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for aldoxime pattern: R-C=N-OH
    aldoxime_pattern = Chem.MolFromSmarts("[#6;R0]=N-O[#1]")
    if mol.HasSubstructMatch(aldoxime_pattern):
        return True, "Contains aldoxime functional group"

    return False, "No aldoxime functional group found"

# Example test
smiles = "CN1C=C(C=C1C=NO)C(=O)C2=CC=CC(=C2)Cl"
print(is_aldoxime(smiles))