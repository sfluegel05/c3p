"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
from rdkit import Chem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    An arenecarbaldehyde is defined as any aldehyde with the carbonyl group attached to an aromatic moiety.

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
    
    # Define an aldehyde pattern with the carbonyl carbon adjacent to an aromatic carbon
    aldehyde_with_aromatic_attachment = Chem.MolFromSmarts("[CX3H1]=[O]")

    # Check if the molecule has an aldehyde group.
    matches = mol.GetSubstructMatches(aldehyde_with_aromatic_attachment)

    for aldehyde_carbon in matches:
        # Find the carbon atom involved in the aldehyde
        aldehyde_atom = mol.GetAtomWithIdx(aldehyde_carbon[0])

        # Check if this carbon atom is attached to an aromatic system
        for neighbor in aldehyde_atom.GetNeighbors():
            if neighbor.GetIsAromatic():
                return True, "Aldehyde group attached to an aromatic moiety"
    
    return False, "Aldehyde group not attached to aromatic moiety"