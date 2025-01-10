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

    # Define SMARTS for potential secosteroid backbone with opened B-ring
    # Include a pattern for a conjugated triene system, common in vitamin D

    # Simplified: triene system would be C=C-C=C-C or similar, should be in any vit-D derivative
    secosteroid_pattern = Chem.MolFromSmarts("[C@H]1C/C=C/C2=CC(CCC2)=C1")

    # Return the classification result based on the SMARTS pattern
    if mol.HasSubstructMatch(secosteroid_pattern):
        return True, "Contains secosteroid structural motif with triene system, suggestive of vitamin D"

    return False, "Does not contain secosteroid structural motif with triene system"

# Example usage
smiles_str = "[C@@H]1(C[C@@H](C/C(=C/C=C/2\CCC[C@]3([C@]2(CC[C@]3([H])[C@](CCCC4(O)CCCC4)([H])C)[H])C)/C1=C)O)O"
result, reason = is_vitamin_D(smiles_str)
print(f"Result: {result}, Reason: {reason}")