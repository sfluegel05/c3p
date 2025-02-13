"""
Classifies: CHEBI:3098 bile acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    A bile acid is a derivative of cholesterol and typically contains a steroid core,
    with various hydroxyl or keto groups, and a side chain ending with a carboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for steroid core pattern (cyclopentanoperhydrophenanthrene ring)
    steroid_core_pattern = Chem.MolFromSmarts("C1C2C3C4C5C(C(C1)C(C23)CC4)CC5")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "No steroid core found"

    # Look for hydroxyl groups attached to the steroid core
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 1:
        return False, "No hydroxyl groups found"

    # Check for keto groups if present and their count
    keto_pattern = Chem.MolFromSmarts("C(=O)")
    keto_matches = mol.GetSubstructMatches(keto_pattern)
    if len(hydroxyl_matches) + len(keto_matches) < 1:
        return False, "Neither sufficient hydroxyl nor keto groups found"

    # Look for carboxylic acid tail
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    return True, "Contains steroid core with oxygen substituents and carboxylic acid side chain"


# Example usage:
smiles_example = "OC(=O)CC[C@H]([C@@]1([C@@]2([C@]([C@]3([C@@]([C@@]4([C@](CC3)(CC=CC4)[H])C)(CC2)[H])[H])(CC1)[H])C)[H])C"
is_bile_acid(smiles_example)