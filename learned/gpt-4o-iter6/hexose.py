"""
Classifies: CHEBI:18133 hexose
"""
from rdkit import Chem

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    Hexoses are defined as six-carbon monosaccharides which in their linear form contain either an aldehyde group at position 1 or a ketone group at position 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check carbon count
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 6:
        return False, f"Incorrect number of carbons for hexose: found {c_count}, expected 6"
    
    # Check for aldehyde group at position 1 (C=O at start of chain)
    aldehyde_pattern = Chem.MolFromSmarts("[CH]=O")
    if mol.HasSubstructMatch(aldehyde_pattern):
        return True, "Contains aldehyde group indicating aldohexose"
    
    # Check for ketone group at position 2 (C-C(=O)-)
    ketone_pattern = Chem.MolFromSmarts("[CH2][C](=O)[CH]")
    if mol.HasSubstructMatch(ketone_pattern):
        return True, "Contains ketone group indicating ketohexose"

    # If neither aldehyde nor ketone pattern matches
    return False, "Does not contain aldehyde or ketone group at required positions for hexose classification"