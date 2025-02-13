"""
Classifies: CHEBI:32876 tertiary amine
"""
from rdkit import Chem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is a tertiary amine based on its SMILES string.
    A tertiary amine has a nitrogen atom bonded to three carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for a tertiary amine (N bonded to 3 carbon atoms)
    tertiary_amine_pattern = Chem.MolFromSmarts("N([C])[C][C]")
    if mol.HasSubstructMatch(tertiary_amine_pattern):
        return True, "Contains a tertiary amine group (N bonded to 3 carbons)"
    else:
        return False, "No tertiary amine group found"

# Example usages
print(is_tertiary_amine("CCN(CC)CC"))  # Expected: True, "Contains a tertiary amine group (N bonded to 3 carbons)"
print(is_tertiary_amine("CCO"))       # Expected: False, "No tertiary amine group found"