"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Define a refined SMARTS pattern for 7-hydroxyisoflavones
    # capturing the benzopyran core specifically and variations
    refined_pattern = Chem.MolFromSmarts("Oc1ccc2c(c1)c(=O)cc(o2)-*")  # Example pattern tweak
    
    # Check if the molecule has the refined 7-hydroxyisoflavone pattern
    if not mol.HasSubstructMatch(refined_pattern):
        return False, "Does not match refined 7-hydroxyisoflavone core structure"
    
    # Check for potential ring closure issues or alternative groups around the core
    # Additional logic can be added to check for variations using further SMARTS or atom checks
    
    # If matched accurately, it is considered a 7-hydroxyisoflavone
    return True, "Matches the refined 7-hydroxyisoflavone core structure"

# Example usage
smiles_examples = [
    "Oc1cc(O)c2c(c1)occ(-c1ccc3OCOc3c1)c2=O",
    "COc1ccc(ccc1O)-c1coc2cc(O)ccc2c1=O",
    # Add more SMILES strings to test
]

for smiles in smiles_examples:
    result, reason = is_7_hydroxyisoflavones(smiles)
    print(f"SMILES: {smiles} -> Is 7-hydroxyisoflavone? {result}. Reason: {reason}")