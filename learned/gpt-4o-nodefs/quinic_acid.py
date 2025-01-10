"""
Classifies: CHEBI:26493 quinic acid
"""
from rdkit import Chem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is a quinic acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Quinic acid core: cyclohexane with multiple oxygen functionalities including a carboxyl group
    # Base pattern with flexibility in hydroxyl groups order and carboxylic acid attachment
    quinic_acid_base_pattern = Chem.MolFromSmarts("O[C@@H]1C[C@@]([C@@H](O)[C@H](O)[C@@H]1O)C(O)=O")
    
    if not mol.HasSubstructMatch(quinic_acid_base_pattern):
        return False, "Mismatch in quinic acid core structure"

    # Additional checks to handle typical quinic acid derivatives could be added here.
    # For current scope, validating the core structure is key.
    
    return True, "Matches quinic acid core structure"