"""
Classifies: CHEBI:35341 steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    A steroid is identified by the cyclopenta[a]phenanthrene carbon skeleton with possible additional structural features.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Basic steroid skeleton pattern (cyclopenta[a]phenanthrene)
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C3CC4CCCC3CCC4C2C1")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid cyclopenta[a]phenanthrene backbone found"
    
    # Check for methyl groups at C-10 and C-13
    methyl_pattern_C10_13 = Chem.MolFromSmarts("[CH3]C1CCC2C3CC4CCCC3CCC4C2C1")
    if not mol.HasSubstructMatch(methyl_pattern_C10_13):
        return False, "Methyl groups at C-10 and/or C-13 are missing"
    
    # Optional: Check for an alkyl group at C-17
    alkyl_pattern_C17 = Chem.MolFromSmarts("[CX4,H]C1CCC2C3CC4CCCC3CCC4C2C1")
    if not mol.HasSubstructMatch(alkyl_pattern_C17):
        return False, "Alkyl group at C-17 is missing"

    # Successful classification
    return True, "Molecule contains a steroid backbone with typical functional groups"

# Examples SMILES test
smiles_examples = [
    "CC(C)C1CCC2C3CC4CCC5=CC(=O)CCC5C4CC3CCC12C",  # Example like cholesterol
    "C[C@@]12[C@](C)C(C)([C@@H]3CC(O)CCC3)(O[C@H](C2)CC1)OC",  # Example like cortisone
    "O[C@@H]1[C@]2(C[C@H](O)[C@](C)(C[C@@H]3[C@H](O)CC[C@]12C3)C)[H]",  # Estrogen derivative
]

for example in smiles_examples:
    print(is_steroid(example))