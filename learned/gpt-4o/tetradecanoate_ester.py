"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
from rdkit import Chem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule is a tetradecanoate ester based on its SMILES string.
    A tetradecanoate ester is obtained by condensation of myristic acid 
    (tetradecanoic acid) with an alcohol or phenol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a tetradecanoate ester, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string to RDKit Mol object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Flexible SMARTS pattern for tetradecanoic acid moiety (considering esters)
    tetradecanoic_acid_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCC(=O)O")

    # Check for tetradecanoic acid moiety as part of an ester
    if not mol.HasSubstructMatch(tetradecanoic_acid_pattern):
        return False, "No tetradecanoic acid moiety found"

    # Check if tetradecanoic acid component is part of an ester linkage
    ester_pattern = Chem.MolFromSmarts("C(=O)O[!H]")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "Tetradecanoic acid moiety not correctly esterified"
    
    return True, "Contains tetradecanoic acid moiety with correct ester linkage"