"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
from rdkit import Chem

def is_aromatic_primary_alcohol(smiles: str):
    """
    Determines if a molecule is an aromatic primary alcohol based on its SMILES string.
    An aromatic primary alcohol is defined as a primary alcohol where the alcoholic hydroxy group 
    is attached to a carbon, which is itself bonded to an aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aromatic primary alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define pattern for an aromatic primary alcohol:
    # [c][CH2][OH] - Aromatic carbon bonded to -CH2- bonded to -OH group
    aromatic_primary_alcohol_pattern = Chem.MolFromSmarts("[c][CH2][OH]")
    
    if aromatic_primary_alcohol_pattern is None:
        return False, "Invalid SMARTS pattern"

    matches = mol.GetSubstructMatches(aromatic_primary_alcohol_pattern)

    if matches:
        for match in matches:
            return True, "Primary alcohol bonded to aromatic ring found"

    return False, "Primary alcohol not bonded to an aromatic ring"