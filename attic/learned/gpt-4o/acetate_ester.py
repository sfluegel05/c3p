"""
Classifies: CHEBI:47622 acetate ester
"""
from rdkit import Chem

def is_acetate_ester(smiles: str):
    """
    Determines if a molecule is an acetate ester based on its SMILES string.
    An acetate ester contains the structural pattern of CH3-C(=O)-O-.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an acetate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define acetate ester pattern (CH3-C(=O)-O-)
    acetate_pattern = Chem.MolFromSmarts("CC(=O)O")
    if acetate_pattern is None:
        return False, "Error in SMARTS pattern creation"

    # Search for the acetate pattern in the molecule
    if mol.HasSubstructMatch(acetate_pattern):
        return True, "Contains acetate ester functional group"
    else:
        return False, "Does not contain acetate ester structural pattern"

# Testing some examples
test_smiles = [
    "CC(=O)OC",    # ethyl acetate
    "CCCCCOC(C)=O",    # pentyl acetate
    "O=C(C)OC",    # methyl acetate
    "CC1=CC(=O)OC1",    # acetoxyethane
]

for smiles in test_smiles:
    result, reason = is_acetate_ester(smiles)
    print(f"SMILES: {smiles} -> Acetate Ester: {result}, Reason: {reason}")