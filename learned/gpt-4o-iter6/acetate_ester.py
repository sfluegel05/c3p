"""
Classifies: CHEBI:47622 acetate ester
"""
from rdkit import Chem

def is_acetate_ester(smiles: str):
    """
    Determines if a molecule is an acetate ester based on its SMILES string.
    An acetate ester contains the esterified acetic acid moiety (CH3COO-).

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

    # Define acetate ester pattern (C(=O)OC)
    acetate_pattern = Chem.MolFromSmarts("C(=O)OC")
    
    # Check if molecule contains the acetate ester pattern
    if mol.HasSubstructMatch(acetate_pattern):
        return True, "Contains acetic acid ester group (CH3COO-)"
    else:
        return False, "Does not contain the acetic acid ester group (CH3COO-)"

# Test examples
smiles_examples = [
    "CC(=O)OC(C(C)C)C(C)C",  # Isopropyl acetate
    "CCCCCCOC(C)=O",          # Hexyl acetate
    "CC1=CC=CC=C1"            # Acetophenone (not an acetate ester)
]

for smiles in smiles_examples:
    result, reason = is_acetate_ester(smiles)
    print(f"SMILES: {smiles}, Is Acetate Ester: {result}, Reason: {reason}")