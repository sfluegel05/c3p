"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
from rdkit import Chem

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    A 2-hydroxy fatty acid contains a carboxylic acid group and a hydroxyl group on the second carbon
    of a hydrocarbon chain.

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

    # Look for the carboxylic acid group pattern as a terminal group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Flexible pattern for hydroxyl group on the second carbon (support for branching/double bonds)
    # The first part `[$(*C(C)O)]` accounts for various chain types (straight/branched with OH)
    hydroxyl_on_second_carbon_pattern = Chem.MolFromSmarts("[$(*C(C)O)]C(=O)O")
    if not mol.HasSubstructMatch(hydroxyl_on_second_carbon_pattern):
        return False, "No hydroxyl group on the second carbon"

    return True, "Contains a carboxylic acid group with a hydroxyl group on the second carbon"