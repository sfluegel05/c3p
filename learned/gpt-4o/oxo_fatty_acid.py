"""
Classifies: CHEBI:59644 oxo fatty acid
"""
from rdkit import Chem

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.
    An oxo fatty acid contains at least one aldehydic or ketonic group in addition to a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an oxo fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group -C(=O)O
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Look for aldehydic group -[CX3H1]=O
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1]=O")
    
    # Look for ketonic group -[CX3](=O)[#6] (must be attached to a carbon)
    ketone_pattern = Chem.MolFromSmarts("[CX3](=O)[#6]")

    # Check for at least one of aldehyde or ketone presence
    has_aldehyde = mol.HasSubstructMatch(aldehyde_pattern)
    has_ketone = mol.HasSubstructMatch(ketone_pattern)

    if not has_aldehyde and not has_ketone:
        return False, "No aldehydic or ketonic group found"

    return True, "Contains both carboxylic acid group and aldehydic or ketonic group, characteristic of an oxo fatty acid"