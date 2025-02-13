"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
"""
Classifies: CHEBI:33762 secondary alpha-hydroxy ketone
"""
from rdkit import Chem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone based on its SMILES string.
    A secondary alpha-hydroxy ketone has a carbonyl group and a hydroxy group linked
    by a carbon bearing one hydrogen and one organyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"

        # Define SMARTS pattern
        secondary_alpha_hydroxy_ketone_smarts = "[CX3H1](=[OX1])[CX3H1]([OX2H1])[CX3]"

        # Check for substructure match
        pattern = Chem.MolFromSmarts(secondary_alpha_hydroxy_ketone_smarts)
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a secondary alpha-hydroxy ketone moiety"
        else:
            return False, "No secondary alpha-hydroxy ketone moiety found"

    except Exception as e:
        return None, f"An error occurred: {str(e)}"