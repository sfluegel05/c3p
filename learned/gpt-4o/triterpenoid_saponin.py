"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid saponin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define triterpenoid backbone pattern (simplified for steroid-like backbone)
    triterpenoid_pattern = Chem.MolFromSmarts("C1CCC2C3CC4CCC(C4)C3CCC2C1")
    if not mol.HasSubstructMatch(triterpenoid_pattern):
        return False, "No triterpenoid backbone found"
    
    # Define glycosidic linkage pattern
    glycoside_pattern = Chem.MolFromSmarts("O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(glycoside_pattern):
        return False, "No glycosidic linkage found"
    
    return True, "Contains triterpenoid structure with glycosidic linkage"

# Test the function with an example
example_smiles = "O[C@H]1[C@H]2C[C@]3(C)[C@H](C[C@@H](O)[C@@H]4[C@H](CC[C@@]34C)[C@@](C)(O)CCC=C(C)C)[C@@]2(C)CC[C@@H]2[C@H]5[C@H](CC[C@@H]1O)[C@@]1(C)[C@@H](O)[C@H](O)[C@H](O)O[C@@H]5C1(C)C2"
result, reason = is_triterpenoid_saponin(example_smiles)
print(result, reason)