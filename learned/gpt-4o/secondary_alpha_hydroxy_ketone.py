"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
from rdkit import Chem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone based on its SMILES string.
    This chemical class involves a hydroxyl group at an alpha position relative to a ketone, 
    where the alpha carbon is a secondary carbon.
    
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
    
    # SMARTS pattern for secondary alpha-hydroxy ketone
    # We'll use:
    # 1. C[OH]: Hydroxyl group bonded to a carbon.
    # 2. C: Secondary nature should be captured by non-terminal, but SMARTS doesn't limit directly; 
    #    hence check hydroxyl presence and alpha connection.
    # 3. C(=O): Alpha position of ketone.
    # Ensuring that the SMARTS match describes a secondary character and ketone relation.
    secondary_alpha_hydroxy_ketone_pattern = Chem.MolFromSmarts('[C@H](O)[CH2][C](=O)')
    
    if not mol.HasSubstructMatch(secondary_alpha_hydroxy_ketone_pattern):
        return False, "No secondary alpha-hydroxy ketone structure found"

    return True, "Contains secondary alpha-hydroxy ketone structure"
    
# Example usage: Test the function with a known secondary alpha-hydroxy ketone.
smiles_example = "O[C@@H]1CC(=CC(=O)[C@H]1O)C(O)=O" # 3-dehydroshikimic acid
result, reason = is_secondary_alpha_hydroxy_ketone(smiles_example)
print(f"Result: {result}, Reason: {reason}")