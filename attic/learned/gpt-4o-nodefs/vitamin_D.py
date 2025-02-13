"""
Classifies: CHEBI:27300 vitamin D
"""
from rdkit import Chem

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D compound based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a vitamin D compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for secosteroid backbone with variable triene configurations
    secosteroid_patterns = [
        # Basic secosteroid backbone with opened B-ring
        Chem.MolFromSmarts("C1[C@@H]2[C@](CCC1=C)([C@@H](CCC(C)C)[H])CCC2"),
        # More flexible pattern that captures triene system post B-ring opening
        Chem.MolFromSmarts("C1CCCC2(C)[C@@H](CCC3=CC=CC=C3)[C@@H](C=C1)C2"),
        # Considerations for variations such as additional functional groups (e.g., hydroxyls)
        Chem.MolFromSmarts("[C@H]([C@H](C)[OH])[C@@H](C)C")
    ]

    # Identify critical structural patterns suggesting vitamin D classification
    found_secosteroid = False
    for pattern in secosteroid_patterns:
        if mol.HasSubstructMatch(pattern):
            found_secosteroid = True
            break
    
    # Determine if it's a candidate for vitamin D based on pattern matching
    if found_secosteroid:
        return True, "Contains structural motifs consistent with vitamin D compounds"
    
    return False, "Does not contain secosteroid structural motifs typical of vitamin D"

# Example usage
smiles_str = "[C@@H]1(C[C@@H](C/C(=C/C=C/2\CCC[C@]3([C@]2(CC[C@]3([H])[C@](CCCC4(O)CCCC4)([H])C)[H])C)/C1=C)O)O"
result, reason = is_vitamin_D(smiles_str)
print(f"Result: {result}, Reason: {reason}")