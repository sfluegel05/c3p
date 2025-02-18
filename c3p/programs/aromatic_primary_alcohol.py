"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
"""
Classifies: CHEBI:134307 aromatic primary alcohol
"""
from rdkit import Chem

def is_aromatic_primary_alcohol(smiles: str):
    """
    Determines if a molecule is an aromatic primary alcohol based on its SMILES string.
    An aromatic primary alcohol has a hydroxyl group (-OH) attached to a primary carbon
    that is directly bonded to an aromatic ring.

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

    # Define SMARTS pattern for primary alcohol attached to aromatic carbon
    # [c] is aromatic carbon, [CH2;D2] is methylene group with exactly two bonds (to OH and aromatic carbon)
    pattern = Chem.MolFromSmarts('[c]-[CH2;D2]-[OH]')
    
    # Check for matches
    matches = mol.GetSubstructMatches(pattern)
    
    if matches:
        return True, "Primary alcohol group attached to aromatic carbon"
    else:
        return False, "No primary alcohol group adjacent to aromatic ring found"