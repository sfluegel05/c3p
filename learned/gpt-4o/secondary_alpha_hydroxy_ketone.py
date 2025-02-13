"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
from rdkit import Chem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone based on its SMILES string.
    A secondary alpha-hydroxy ketone has an alpha carbon bonded to one hydroxyl group, 
    one hydrogen atom and two carbons (making it a secondary carbon) that is adjacent 
    to a carbonyl (C=O) group.
    
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
    
    # SMARTS for secondary alpha-hydroxy ketone
    # `[C;H1&R0]`: Secondary carbon (H1 means single hydrogen, indicating secondary, and R0 to ensure no cyclic reference)
    # `[OH]`: directly attached hydroxyl group
    # `[C]=[O]`: carbonyl group adjacent to the alpha position
    secondary_alpha_hydroxy_ketone_pattern = Chem.MolFromSmarts('[C;H1&R0][C;H1]([OH])[R0][C]=[O]')
    if not mol.HasSubstructMatch(secondary_alpha_hydroxy_ketone_pattern):
        return False, "No secondary alpha-hydroxy ketone structure found"

    return True, "Contains secondary alpha-hydroxy ketone structure"

# Example usage: Test the function with a known secondary alpha-hydroxy ketone.
smiles_example = "OC[C@@H]1CC=C(O)C(=O)[C@H]1O" # 5D-(5/6)-2,6-dihydroxy-5-(hydroxymethyl)cyclohex-2-en-1-one
result, reason = is_secondary_alpha_hydroxy_ketone(smiles_example)
print(f"Result: {result}, Reason: {reason}")