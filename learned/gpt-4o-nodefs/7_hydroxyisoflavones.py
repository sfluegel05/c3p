"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
from rdkit import Chem

def is_7_hydroxyisoflavones(smiles: str):
    """
    Determines if a molecule is a 7-hydroxyisoflavone based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 7-hydroxyisoflavone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for 7-hydroxyisoflavones
    # capturing the isoflavone core and ensuring a hydroxyl group at the 7-position
    refined_pattern = Chem.MolFromSmarts("Oc1ccc2c(c1)occ(-c1ccc(O)cc1)c2=O")
    
    # Check if the molecule has the refined 7-hydroxyisoflavone pattern
    if not mol.HasSubstructMatch(refined_pattern):
        return False, "Does not match refined 7-hydroxyisoflavone core structure"
    
    # Check for any alternative functional groups or additional variability around the isoflavone core
    # Logic can be enhanced by adding additional checks if specific derivatizations are expected
    
    return True, "Matches the refined 7-hydroxyisoflavone core structure"

# Example usage
smiles_examples = [
    "Oc1cc(O)c2c(c1)occ(-c1ccc3OCOc3c1)c2=O",  # Example 7-hydroxyisoflavone
    "COc1ccc(ccc1O)-c1coc2cc(O)ccc2c1=O",  # Example
    # Add more SMILES strings for testing
]

for smiles in smiles_examples:
    result, reason = is_7_hydroxyisoflavones(smiles)
    print(f"SMILES: {smiles} -> Is 7-hydroxyisoflavone? {result}. Reason: {reason}")