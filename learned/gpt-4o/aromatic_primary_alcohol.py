"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
from rdkit import Chem

def is_aromatic_primary_alcohol(smiles: str):
    """
    Determines if a molecule is an aromatic primary alcohol based on its SMILES string.
    An aromatic primary alcohol has a hydroxy group attached to a carbon which is itself bonded to an aromatic ring.

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

    # SMARTS pattern for a primary alcohol where the OH is bonded to sp3 carbon
    # and that sp3 carbon is bonded to an aromatic carbon
    pattern = Chem.MolFromSmarts("[CX4;$([CH2]O)]-[OH]-[a]")

    # Check for the presence of the aromatic primary alcohol pattern
    if mol.HasSubstructMatch(pattern):
        return True, "Contains hydroxy group attached to a carbon connected to an aromatic ring"
    
    return False, "Does not contain aromatic primary alcohol structure"

# The function can now be used to classify SMILES strings for aromatic primary alcohols.