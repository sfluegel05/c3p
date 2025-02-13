"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
from rdkit import Chem

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.
    An alpha-hydroxy ketone contains a ketone group with a hydroxy group on the alpha-carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for combined pattern of ketone with alpha-hydroxy group
    alpha_hydroxy_ketone_pattern = Chem.MolFromSmarts("[C;H1,H2](O)[C;H2,H3](=O)[C]")
    if not mol.HasSubstructMatch(alpha_hydroxy_ketone_pattern):
        return False, "No alpha-hydroxy ketone structure found"

    return True, "Contains a ketone group with a hydroxy group on the alpha-carbon"