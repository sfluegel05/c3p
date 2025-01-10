"""
Classifies: CHEBI:20706 6-aminopurines
"""
"""
Classifies: 6-aminopurines
Definition: Any compound having 6-aminopurine (adenine) as part of its structure
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_6_aminopurines(smiles: str):
    """
    Determines if a molecule contains a 6-aminopurine (adenine) moiety based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains 6-aminopurine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for 6-aminopurine (adenine)
    # This pattern represents the fused bicyclic system with amino group at position 6
    # The pattern matches:
    # - A 5-membered ring fused to a 6-membered ring (purine scaffold)
    # - An amino group (-NH2) at position 6
    # - The correct number and position of nitrogens in the rings
    adenine_pattern = Chem.MolFromSmarts('c1[nH]c2c(n1)c(N)nc[nH]2')
    
    # Alternative pattern that also matches N9-substituted adenines (as in nucleotides)
    adenine_pattern_2 = Chem.MolFromSmarts('c1nc2c(n1)c(N)ncn2')
    
    # Check for matches
    has_adenine_1 = mol.HasSubstructMatch(adenine_pattern)
    has_adenine_2 = mol.HasSubstructMatch(adenine_pattern_2)
    
    if has_adenine_1 or has_adenine_2:
        return True, "Contains 6-aminopurine (adenine) moiety"
    
    return False, "Does not contain 6-aminopurine structure"