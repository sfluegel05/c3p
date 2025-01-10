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

    # Define a more general SMARTS pattern for a secosteroid backbone
    # with opened B-ring and variations of conjugated triene
    # Note: This pattern may require further adjustment based on structural diversity
    # of Vitamin D derivatives. Adjusted to include characteristic open B-ring and triene systems.
    secosteroid_patterns = [
        Chem.MolFromSmarts("[C@H]1C(=C)[CCC@@]2CC=C1C3=CC=C[C@]3(C)C4([C@@H]2CCC4)C"), # typical secosteroid
        Chem.MolFromSmarts("C[C@H]1[C@H](C)[C@@H](C/C=C\\C2=C/[C@H](CC[C@@](C3)(C(=CCC3)=C1)C\\C)CC2)=CC") # triene with opened B-ring
        # More patterns as needed based on common forms seen in vitamin D
    ]

    # Check each pattern
    for pattern in secosteroid_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains secosteroid structural motif, suggestive of vitamin D"

    return False, "Does not contain secosteroid structural motif"

# Example usage
smiles_str = "[C@@H]1(C[C@@H](C/C(=C/C=C/2\CCC[C@]3([C@]2(CC[C@]3([H])[C@](CCCC4(O)CCCC4)([H])C)[H])C)/C1=C)O)O"
result, reason = is_vitamin_D(smiles_str)
print(f"Result: {result}, Reason: {reason}")