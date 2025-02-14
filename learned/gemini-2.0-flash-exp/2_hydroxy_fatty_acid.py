"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    A 2-hydroxy fatty acid has a carboxylic acid group and a hydroxyl group
    on the carbon alpha to the carboxyl carbon.

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

    # Check for the presence of a carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("C(=O)[O,OH]")
    if not mol.HasSubstructMatch(acid_pattern):
         return False, "No carboxylic acid group found"

    # Check for 2-hydroxy group (hydroxy on the alpha carbon)
    hydroxy_pattern = Chem.MolFromSmarts("[CX4,CX3]([OH])C(=O)[OH,O]")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No hydroxyl group at the alpha position found"

    return True, "Contains a carboxylic acid group with a hydroxyl group on the alpha carbon"