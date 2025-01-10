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
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Isoflavone core pattern with a hydroxy group at the 7th position
    # Note: This is a very simplified pattern to capture the core structure
    hydroxyisoflavone_pattern = Chem.MolFromSmarts("Oc1ccc2oc(=O)ccc2c1")
    
    # Check for the 7-hydroxyisoflavone pattern
    if not mol.HasSubstructMatch(hydroxyisoflavone_pattern):
        return False, "Does not match 7-hydroxyisoflavone core structure"
    
    # If matched, then it's a 7-hydroxyisoflavone
    return True, "Matches the 7-hydroxyisoflavone core structure"

# Example usage
smiles_examples = [
    "Oc1cc(O)c2c(c1)occ(-c1ccc3OCOc3c1)c2=O",
    "COc1ccc(ccc1O)-c1coc2cc(O)ccc2c1=O",
    # Add more SMILES strings to test
]

for smiles in smiles_examples:
    result, reason = is_7_hydroxyisoflavones(smiles)
    print(f"SMILES: {smiles} -> Is 7-hydroxyisoflavone? {result}. Reason: {reason}")