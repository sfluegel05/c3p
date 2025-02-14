"""
Classifies: CHEBI:134363 tertiary amine oxide
"""
from rdkit import Chem

def is_tertiary_amine_oxide(smiles: str):
    """
    Determines if a molecule is a tertiary amine oxide based on its SMILES string.
    
    A tertiary amine oxide is characterized by a positively charged nitrogen atom 
    bonded to three organic groups, with one of its lone pairs forming a bond with an oxygen atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine oxide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for a tertiary amine oxide pattern
    # This pattern looks for an N-oxide where nitrogen has three organic groups attached
    amine_oxide_pattern = Chem.MolFromSmarts("[N+]([C,R])([C,R])([C,R])[O-]")
    
    if mol.HasSubstructMatch(amine_oxide_pattern):
        return True, "Contains tertiary amine oxide structure"
    else:
        return False, "Does not contain tertiary amine oxide structure"

# Example usage
smiles = "C[N+](C)([O-])C"  # Trimethylamine N-oxide
result, reason = is_tertiary_amine_oxide(smiles)
print(f"Result: {result}, Reason: {reason}")