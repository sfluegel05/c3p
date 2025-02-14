"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
from rdkit import Chem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    An arenecarbaldehyde is defined as an aldehyde in which the carbonyl group is attached to an aromatic moiety.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an arenecarbaldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define aldehyde group pattern
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1]=O")
    
    # Check if the molecule contains an aldehyde group
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)

    # Define pattern for aldehyde attached to an aromatic carbon
    aromatic_aldehyde_pattern = Chem.MolFromSmarts("c[CX3H1]=O")

    # Check for aldehyde attached to aromatic carbon
    for match in aldehyde_matches:
        aldehyde_carbon = match[0]
        # Check if the aldehyde carbon is directly connected to an aromatic carbon
        for neighbor in mol.GetAtomWithIdx(aldehyde_carbon).GetNeighbors():
            if neighbor.GetIsAromatic():
                return True, "Aromatic ring with attached aldehyde group found"

    return False, "No aromatic ring directly bound to an aldehyde group found"