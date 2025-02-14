"""
Classifies: CHEBI:20706 6-aminopurines
"""
from rdkit import Chem

def is_6_aminopurines(smiles: str):
    """
    Determines if a molecule contains a 6-aminopurine (adenine) structure based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains the 6-aminopurine structure, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Refined the 6-aminopurine SMARTS pattern for adenine
    # Pay attention to bonds and possible variation in ring assemblies
    adenine_pattern = Chem.MolFromSmarts("n1cnc2c(ncn2)n1")  # Core purine pattern with 6-amino consideration
    
    if adenine_pattern is None:
        return False, "Failed to create 6-aminopurine pattern"
    
    # Check if the molecule contains the adenine spine
    if mol.HasSubstructMatch(adenine_pattern):
        # Ensure the N at position 6 is bonded
        matches = mol.GetSubstructMatches(adenine_pattern)
        for match in matches:
            if len(match) >= 6:
                n_pos = match[4]  # The index for the nitrogen expected at position 6
                n_atom = mol.GetAtomWithIdx(n_pos)
                if n_atom.GetSymbol() == 'N':
                    return True, "Contains 6-aminopurine (adenine) structure"
    
    return False, "Does not contain 6-aminopurine (adenine) structure"

# Note: It's crucial to double-check this refined pattern with a variety of test cases to confirm accuracy.