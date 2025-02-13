"""
Classifies: CHEBI:50523 butenolide
"""
from rdkit import Chem

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide based on its SMILES string.
    A butenolide consists of a 5-membered lactone ring (furanone structure), often with a ketone group or 
    related functional group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butenolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for butenolide structure
    # More flexible pattern to catch variations in the butenolide structure
    butenolide_pattern_1 = Chem.MolFromSmarts("C1=COC(=O)C1")  # Basic furanone pattern
    butenolide_pattern_2 = Chem.MolFromSmarts("C1=CC(=O)OC1")  # Variation on furanone positioning

    # Check if molecule matches any of the butenolide substructures
    if mol.HasSubstructMatch(butenolide_pattern_1):
        return True, "Contains a basic furanone (butenolide) structure"
    elif mol.HasSubstructMatch(butenolide_pattern_2):
        return True, "Contains a variation of furanone (butenolide) structure"
    
    return False, "Does not contain a recognized butenolide structure"

# Test the function with examples
examples = [
    "CC1(C)OC(=O)C(OCC2CC2)=C1c1ccc(cc1)S(C)(=O)=O",  # firocoxib
    "O=C1OC(=C)C(=C1C(=O)CCCC=CCCCCCC)O",            # Agglomerin B
]

for example in examples:
    result, reason = is_butenolide(example)
    print(f"SMILES: {example} => Result: {result}, Reason: {reason}")