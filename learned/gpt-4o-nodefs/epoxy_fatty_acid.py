"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
"""
Classifies: epoxy fatty acid
"""
from rdkit import Chem

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.
    An epoxy fatty acid features an epoxide group and a carboxylic acid group attached to a hydrocarbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an epoxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify epoxide group
    epoxide_pattern = Chem.MolFromSmarts("C1OC1")
    if not mol.HasSubstructMatch(epoxide_pattern):
        return False, "No epoxide group found"
    
    # Identify carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Identify long hydrocarbon chain
    chain_pattern = Chem.MolFromSmarts("CCCCCC")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No long hydrocarbon chain found"
        
    return True, "Contains epoxide group and carboxylic acid along with a hydrocarbon chain"