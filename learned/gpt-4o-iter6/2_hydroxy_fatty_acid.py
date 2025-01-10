"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
from rdkit import Chem

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    A 2-hydroxy fatty acid has a hydroxy group in the alpha- or 2-position relative to the carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 2-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define pattern for a 2-hydroxy fatty acid with potential stereo-centers
    # [C@@H]- alcohol for indicating R or S stereocenter; |!@ (no stereochemistry) accounts for linear connectivity
    # [-][O]C(O)=O for the carboxylic acid pattern
    pattern = Chem.MolFromSmarts("[CX4,CX3](O)[CX3](=O)O")

    if mol.HasSubstructMatch(pattern):
        return True, "Contains 2-hydroxy group with carboxylic acid at the 2-position"
    
    return False, "No 2-hydroxy group found at the alpha position"

# Example usage
print(is_2_hydroxy_fatty_acid("CCCCCCC[C@@H](O)C(O)=O"))  # (R)-2-hydroxynonanoic acid