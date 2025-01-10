"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid based on its SMILES string.
    A hydroxy fatty acid is defined as a fatty acid carrying one or more hydroxy substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group pattern
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O,H]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Look for hydroxy group (-OH)
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No hydroxy groups found"
    
    return True, "Contains carboxylic acid group with one or more hydroxy substituents"

__metadata__ = {
    'chemical_class': {
        'name': 'hydroxy fatty acid',
        'definition': 'Any fatty acid carrying one or more hydroxy substituents.'
    }
}