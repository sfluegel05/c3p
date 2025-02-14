"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
"""
Classifies: CHEBI:136568 2-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    A 2-hydroxy fatty acid is any fatty acid with a hydroxy functional group in the alpha- or 2-position.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a 2-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for 2-hydroxy fatty acid
    # Carboxylic acid group connected to a carbon with hydroxy group (alpha carbon)
    pattern = Chem.MolFromSmarts('[CX3](=O)[OX1H1][CH](O)[C,c]')
    if pattern is None:
        return False, "Invalid SMARTS pattern"

    # Check for substructure match
    if mol.HasSubstructMatch(pattern):
        return True, "Matches 2-hydroxy fatty acid structure"
    else:
        return False, "Does not match 2-hydroxy fatty acid structure"