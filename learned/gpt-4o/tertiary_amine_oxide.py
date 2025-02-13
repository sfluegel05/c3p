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
    
    # Define tertiary amine oxide pattern [Nitrogen with three organic groups, one O-]
    pattern = Chem.MolFromSmarts("[NX4+](~[O-])(~[C])(~[C])[C]")
    if mol.HasSubstructMatch(pattern):
        return True, "Contains tertiary amine oxide structure with correct arrangement"
        
    return False, "Does not contain the correct tertiary amine oxide structure"

# Example usage
test_smiles = "C[N+](C)([O-])C"  # Trimethylamine N-oxide
result, reason = is_tertiary_amine_oxide(test_smiles)
print(f"Result: {result}, Reason: {reason}")