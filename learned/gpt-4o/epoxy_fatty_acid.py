"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
from rdkit import Chem

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.
    An epoxy fatty acid contains an epoxide ring as part of its structure,
    along with a long aliphatic chain typically ending in a carboxyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an epoxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of an epoxide group (cyclic ether C1OC1)
    epoxide_pattern = Chem.MolFromSmarts("C1OC1")
    if not mol.HasSubstructMatch(epoxide_pattern):
        return False, "No epoxide group found"
    
    # Check for the carboxylic acid group at the end of a long aliphatic chain
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")  # Simplified pattern for carboxylic acid
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for the presence of the long aliphatic chain
    aliphatic_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCC")  # At least 10 carbon chain pattern
    if not mol.HasSubstructMatch(aliphatic_chain_pattern):
        return False, "Aliphatic chain too short to be a fatty acid"
    
    return True, "Contains epoxide ring as part of a long aliphatic chain ending in a carboxyl group"