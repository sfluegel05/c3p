"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
from rdkit import Chem

def is_aromatic_primary_alcohol(smiles: str):
    """
    Determines if a molecule is an aromatic primary alcohol based on its SMILES string.
    An aromatic primary alcohol is defined as a primary alcohol where the alcoholic hydroxy group 
    is attached to a carbon which is itself bonded to an aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic primary alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define primary alcohol attached to aromatic pattern:
    # It is [CH2][OH] where [CH2] connects to an aromatic ring [c].
    aromatic_primary_alcohol_pattern = Chem.MolFromSmarts("[c][CH2][OH]")

    # Find substructures matching the aromatic primary alcohol pattern
    if aromatic_primary_alcohol_pattern is None:
        return False, "Invalid SMARTS pattern"

    matches = mol.GetSubstructMatches(aromatic_primary_alcohol_pattern)
    
    if matches:
        return True, "Primary alcohol bonded to aromatic ring found"

    return False, "Primary alcohol not bonded to an aromatic ring"

# Example usage
# smile = "CC1=C(CO)C=CC=C1"  # 2-methylbenzyl alcohol
# print(is_aromatic_primary_alcohol(smile))