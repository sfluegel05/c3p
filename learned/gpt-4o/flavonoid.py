"""
Classifies: CHEBI:47916 flavonoid
"""
from rdkit import Chem

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    A flavonoid is defined broadly as a 1-benzopyran with an aryl substituent at position 2,
    allowing for various decorations on the flavonoid core structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern with more flexibility around the flavonoid core structure and an aryl group at position 2
    flavonoid_smarts = "[O]1c2ccc(O)cc2Oc3ccccc3-c1"  # captures a typical flavonoid core with flexibility for common substituents
    
    flavonoid_pattern = Chem.MolFromSmarts(flavonoid_smarts)

    # Check for flavonoid structure
    if not mol.HasSubstructMatch(flavonoid_pattern):
        return False, "No flavonoid structure detected"

    return True, "Molecule contains the core 1-benzopyran structure with aryl group typical of flavonoids"

# Test with example SMILES
sample_smiles = [
    "O1C2=C(C(OC)=C3C(OC=C3)=C2)C(=O)C=C1C4=CC=CC=C4",  # Pinnatin
    "O1C2=C(CC=C(C)C)C(O)=CC(O)=C2C(=O)C=C1C3=CC(OC)=C(O)C(OC)=C3",  # Baohuosu
    "O1C([C@@H](O)[C@H](O)C(O)[C@@H]1OC=2C(=[O+]C=3C(C2)=C(O)C=C(O)C3)C4=CC(O)=C(O)C(O)=C4)CO[C@@H]5OC([C@H](O)C(O)[C@@H]5OC(=O)C)C",  # Delphinidin 3-[6-(2-acetylrhamnosyl)glucoside]
]

# Validate with the examples
results = [(smiles, is_flavonoid(smiles)) for smiles in sample_smiles]
for smiles, (result, reason) in results:
    print(f"SMILES: {smiles}\nResult: {result}\nReason: {reason}\n")