"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
from rdkit import Chem

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.
    An alpha-hydroxy ketone has a hydroxyl group directly adjacent to a ketone group (C=O),
    accounting for both acyclic and cyclic structures.

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
    
    # Improve the alpha-hydroxy ketone pattern by considering both acyclic and cyclic structures
    alpha_hydroxy_ketone_pattern_general = Chem.MolFromSmarts("[$(OC(=O));!$(*=*])")  # simple alpha-hydroxy-ketone motif
    alpha_hydroxy_ketone_pattern_cyclic = Chem.MolFromSmarts("[$(O)C(=O)]")  # general cyclic

    # Check if the molecule has the alpha-hydroxy ketone substructure pattern
    if mol.HasSubstructMatch(alpha_hydroxy_ketone_pattern_general) or mol.HasSubstructMatch(alpha_hydroxy_ketone_pattern_cyclic):
        return True, "Contains an alpha-hydroxy ketone group"
    else:
        return False, "Does not contain an alpha-hydroxy ketone group"

# Example usage:
# smiles = "CC(=O)CO" # hydroxyacetone
# is_alpha_hydroxy_ketone(smiles)