"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
"""
Classifies: aliphatic alcohol (CHEBI:50584)
"""
from rdkit import Chem

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.
    An aliphatic alcohol has at least one hydroxyl group (-OH) attached to an aliphatic carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic alcohol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Sanitize to detect aromaticity
    Chem.SanitizeMol(mol)
    
    # Find all hydroxyl groups (-OH) attached to a carbon
    # SMARTS pattern: oxygen with one H connected to any carbon
    pattern = Chem.MolFromSmarts('[O;H1][C]')
    matches = mol.GetSubstructMatches(pattern)
    
    if not matches:
        return False, "No hydroxyl group found"
    
    # Check if any hydroxyl is attached to a non-aromatic (aliphatic) carbon
    for match in matches:
        oxygen_idx, carbon_idx = match
        carbon = mol.GetAtomWithIdx(carbon_idx)
        if not carbon.GetIsAromatic():
            return True, "Hydroxyl group attached to aliphatic carbon"
    
    # All hydroxyls are on aromatic carbons
    return False, "All hydroxyl groups are on aromatic carbons"