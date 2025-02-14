"""
Classifies: CHEBI:35341 steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    A steroid is identified by the cyclopenta[a]phenanthrene carbon skeleton with flexible pattern matching
    for typical structural variations and functional groups.
    
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

    # Generalized steroid backbone pattern - cyclopenta[a]phenanthrene with flexible saturation
    steroid_pattern = Chem.MolFromSmarts("C12C3C4C1C5C3C(C4)CCC5C2")
    if mol.HasSubstructMatch(steroid_pattern):
        # Optional: Check for typical methyl group presence
        methyl_checker = Chem.MolFromSmarts("[C@@H]([CH3])[C@H]1C[C@@H]2CC[C@]3C1[C@H]2CC[C@@]3C")
        methyl_matches = mol.GetSubstructMatches(methyl_checker)
        if methyl_matches:
            return True, "Molecule contains a steroid backbone with typical methyl groups."

        return True, "Molecule contains a steroid backbone."

    return False, "No detectable steroid backbone structure found."

# Examples SMILES test
smiles_examples = [
    "C[C@]12CC[C@H]3[C@H]([C@H]1CC[C@]2(C#C)O)CCC4=C3C=CC(=C4)OC",  # Expected True
    "CC(C)=CCC[C@](C)(O)[C@H]1CC[C@]2(C)[C@@H]1[C@H](O)C[C@@H]1[C@@]3(C)CC[C@H](O)C(C)(C)[C@@H]3CC[C@@]21C",  # True
]

for example in smiles_examples:
    print(is_steroid(example))