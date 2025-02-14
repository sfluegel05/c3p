"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.
    A 3-hydroxy fatty acid has a hydroxyl group at the beta- or 3-position relative to the carboxyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify the pattern: a carboxylic acid attached to a carbon chain with an OH on the 3rd carbon
    # The SMARTS pattern captures [C]-[C](-O)[C(=O)O], where the first [C] is the first
    # carbon chain position, the second [C] is the 2nd position with a hydroxyl group, and [C(=O)O] is the carboxyl group.
    hydroxy_fatty_acid_pattern = Chem.MolFromSmarts("CCO[C@H][C](=O)O")
    
    # Check for presence of the 3-hydroxy fatty acid pattern
    if mol.HasSubstructMatch(hydroxy_fatty_acid_pattern):
        return True, "Contains hydroxyl group at the 3-position and carboxylic acid group"

    # Adding handling for potential stereoisomer patterns, i.e., the presence of @ or @@ in SMILES
    stereo_hydroxy_fatty_acid_pattern = Chem.MolFromSmarts("C[C@H](O)[C](=O)O")
    if mol.HasSubstructMatch(stereo_hydroxy_fatty_acid_pattern):
        return True, "Contains stereochemical hydroxyl group pattern at the 3-position and carboxylic acid group"

    return False, "Hydroxyl group not found at the 3-position or missing carboxyl group"