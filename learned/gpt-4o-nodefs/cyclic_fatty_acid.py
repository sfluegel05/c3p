"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid based on its SMILES string.
    Cyclic fatty acids generally contain a cycle within a long aliphatic chain and a carboxylic acid group.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclic fatty acid, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of a carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for the presence of at least one ring
    if not rdmolops.GetSSSR(mol):
        return False, "No cyclic structure detected"
    
    # Check for sufficient length of carbon chain
    carbon_chain_length = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_chain_length < 10:
        return False, f"Carbon chain too short: {carbon_chain_length} carbons"

    # Check for common cyclic fatty acid features
    # For example, look for epoxide substructures if applicable
    epoxide_pattern = Chem.MolFromSmarts("[OX2][CX4][OX2]")
    has_epoxide = mol.HasSubstructMatch(epoxide_pattern)

    # Additional checks can be implemented as needed
    # ...

    return True, "Contains a cyclic structure, carboxylic acid group, and sufficient carbon chain length."