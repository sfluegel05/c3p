"""
Classifies: CHEBI:50523 butenolide
"""
from rdkit import Chem

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide based on its SMILES string.
    A butenolide consists of a 5-membered lactone ring, or furanone structure, with a ketone group.

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
    butenolide_pattern = Chem.MolFromSmarts("O=C1C=CC(=O)O1")
    
    # Check if molecule has a butenolide substructure
    if mol.HasSubstructMatch(butenolide_pattern):
        return True, "Contains a butenolide (furanone) structure"
    
    return False, "Does not contain a butenolide structure"

# Test the function with examples
examples = [
    "CC1(C)OC(=O)C(OCC2CC2)=C1c1ccc(cc1)S(C)(=O)=O",  # firocoxib
    "O=C1OC(=C)C(=C1C(=O)CCCC=CCCCCCC)O",            # Agglomerin B
]

for example in examples:
    result, reason = is_butenolide(example)
    print(f"SMILES: {example} => Result: {result}, Reason: {reason}")