"""
Classifies: CHEBI:35342 17alpha-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_17alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17alpha-hydroxy steroid based on its SMILES string.
    A 17alpha-hydroxy steroid has a hydroxyl group at the 17alpha position on the steroid structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for steroid backbone (cyclopentanoperhydrophenanthrene core)
    steroid_pattern = Chem.MolFromSmarts('C1CCC2C3CCC4C(C3C2C1)CC4')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Define SMARTS pattern for the 17alpha-hydroxy group
    hydroxyl_pattern = Chem.MolFromSmarts('[C@H](O)[C&R3,R4]')
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No 17alpha-hydroxyl group found"
    
    # Additional logic can be added here to check stereochemistry more accurately if needed

    return True, "Contains steroid backbone with 17alpha-hydroxyl group"

# Example usage
smiles_example = '[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])[C@@H](O)C[C@@]1(C)[C@@]2([H])CC[C@]1(O)C(=O)CO'
print(is_17alpha_hydroxy_steroid(smiles_example))