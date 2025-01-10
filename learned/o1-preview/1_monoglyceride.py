"""
Classifies: CHEBI:35759 1-monoglyceride
"""
from rdkit import Chem

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 1-monoglyceride based on its SMILES string.
    A 1-monoglyceride is a monoglyceride in which the acyl substituent is located at position 1.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a 1-monoglyceride, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the monoglyceride pattern with ester at position 1
    # Position 1: esterified hydroxyl
    # Position 2: secondary alcohol
    # Position 3: primary alcohol
    monoglyceride_pattern = Chem.MolFromSmarts("C(=O)O[CH2][CH](O)[CH2]O")
    
    # Check for substructure match without considering chirality
    if mol.HasSubstructMatch(monoglyceride_pattern, useChirality=False):
        return True, "Molecule matches 1-monoglyceride structure with acyl group at position 1 and free hydroxyls at positions 2 and 3"

    return False, "Molecule does not match 1-monoglyceride structure"