"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
from rdkit import Chem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone based on its SMILES string.
    A secondary alpha-hydroxy ketone has a carbonyl group and hydroxy group linked by a carbon bearing one hydrogen and one organyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for secondary alpha-hydroxy ketone
    secondary_alpha_hydroxy_ketone_pattern = Chem.MolFromSmarts("[CX4;H1][CX3](=O)[OX2H]")
    
    # Check for the presence of the pattern
    if mol.HasSubstructMatch(secondary_alpha_hydroxy_ketone_pattern):
        return True, "Molecule contains the characteristic structure of a secondary alpha-hydroxy ketone"
    
    return False, "No characteristic secondary alpha-hydroxy ketone structure found"

# Test examples with known secondary alpha-hydroxy ketones
test_smiles = "O=C1C=C[C@H]([C@]1(O)C(=O)CCCCCCCCCCCC)O"  # Hygrophorone C12
is_secondary_alpha_hydroxy_ketone(test_smiles)