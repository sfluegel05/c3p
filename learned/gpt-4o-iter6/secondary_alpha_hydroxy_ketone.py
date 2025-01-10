"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
from rdkit import Chem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone based on its SMILES string.
    A secondary alpha-hydroxy ketone has an -OH group adjacent to a C=O group
    where the central carbon is secondary (has one hydrogen and one other organic group).

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

    # Look for a secondary alpha-hydroxy ketone pattern, C with OH and adjacent C=O, C must be secondary
    # C=O group [C]=[O] adjacent to a [C][OH] where C has one H and one organic (R) group, no other specifics needed
    pattern = Chem.MolFromSmarts("[#6;R1;$([CD3](=O)[C,O]);$([CD3]([OX2H])([C][H]))]")
    
    # Check for match of the pattern
    if mol.HasSubstructMatch(pattern):
        return True, "Contains a secondary alpha-hydroxy ketone functional group"
        
    return False, "Does not contain a secondary alpha-hydroxy ketone functional group"