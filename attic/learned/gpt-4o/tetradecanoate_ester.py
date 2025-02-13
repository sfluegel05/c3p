"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
from rdkit import Chem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule is a tetradecanoate ester
    based on its SMILES string.
    A tetradecanoate ester is obtained by condensation of
    myristic acid (tetradecanoic acid) with an alcohol or phenol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetradecanoate ester, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string to RDKit Mol object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the tetradecanoic acid moiety pattern using SMILES
    tetradecanoic_acid_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCC(=O)O")
    # Define ester linkage pattern
    ester_linkage_pattern = Chem.MolFromSmarts("C(=O)O[C,N]")

    # Check for tetradecanoic acid moiety
    if not mol.HasSubstructMatch(tetradecanoic_acid_pattern):
        return False, "No tetradecanoic acid moiety found"

    # Check for ester linkage
    if not mol.HasSubstructMatch(ester_linkage_pattern):
        return False, "No ester linkage found"

    return True, "Contains tetradecanoic acid moiety with ester linkage"