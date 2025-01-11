"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
from rdkit import Chem

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid based on its SMILES string.
    An unsaturated fatty acid contains at least one C=C or C#C bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an unsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least one C=C or C#C bond
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    triple_bond_pattern = Chem.MolFromSmarts("C#C")
    
    if not mol.HasSubstructMatch(double_bond_pattern) and not mol.HasSubstructMatch(triple_bond_pattern):
        return False, "No C=C or C#C bonds found"

    # Check for carboxylic acid group to verify it's a fatty acid
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found, not a fatty acid"
    
    return True, "Contains at least one C=C or C#C bond and a carboxylic acid group"