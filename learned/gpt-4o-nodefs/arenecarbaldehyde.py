"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
from rdkit import Chem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    An arenecarbaldehyde is characterized by an aldehyde group directly attached to an aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an arenecarbaldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for aromatic ring (any kind) with aldehyde group directly attached
    arenecarbaldehyde_pattern = Chem.MolFromSmarts("[c][C](=O)")

    # Verify the molecule contains the necessary substructure
    if mol.HasSubstructMatch(arenecarbaldehyde_pattern):
        # Additional verification: Confirm that the entire fragment behaves as desired in aromatic context
        submatches = mol.GetSubstructMatches(arenecarbaldehyde_pattern)
        for match in submatches:
            aldehyde_carbon = match[1]  # Index of the aldehyde carbon based on the SMARTS pattern
            neighbors = mol.GetAtomWithIdx(aldehyde_carbon).GetNeighbors()
            for neighbor in neighbors:
                if neighbor.GetIsAromatic():
                    return True, "Contains aromatic ring with an aldehyde group directly attached"
        return False, "Pattern found but not in correct aromatic context"
    else:
        return False, "No arenecarbaldehyde structure found"