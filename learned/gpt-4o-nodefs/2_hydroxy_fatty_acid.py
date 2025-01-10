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
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxylic_acid_matches:
        return False, "No carboxylic acid group found"

    # Ensure the hydroxyl is on the second carbon from the carboxylic acid
    # To address branching and other configurations flexibility is added in substructure matching
    hydroxy_on_second_carbon_pattern = Chem.MolFromSmarts("[C;H1,H2,H3][C;H1,H2](O)C(=O)O")
    if not mol.HasSubstructMatch(hydroxy_on_second_carbon_pattern):
        return False, "No hydroxyl group on the second carbon from the carboxylic acid"

    # Additional checks for chain length and unsaturation depend on further specifics not covered here
    return True, "Contains carboxylic acid and hydroxyl group on the second carbon as expected in 2-hydroxy fatty acids"