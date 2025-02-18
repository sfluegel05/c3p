"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
"""
Classifies: arenecarbaldehyde
"""

from rdkit import Chem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    An arenecarbaldehyde is any aldehyde in which the carbonyl group is attached to an aromatic moiety.
    
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
    
    # Define the substructure pattern for aldehyde group
    aldehyde_pattern = Chem.MolFromSmarts('[CH]=O')
    if aldehyde_pattern is None:
        return False, "Failed to create aldehyde pattern"
    
    # Find all aldehyde groups in the molecule
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    n_aldehydes = len(aldehyde_matches)
    
    if n_aldehydes == 0:
        return False, "No aldehyde groups found in the molecule"
    
    # Define the substructure pattern for aldehyde attached to aromatic ring
    aromatic_aldehyde_pattern = Chem.MolFromSmarts('[a][CH]=O')
    if aromatic_aldehyde_pattern is None:
        return False, "Failed to create aromatic aldehyde pattern"
        
    # Find all aldehyde groups attached to aromatic moieties
    aromatic_aldehyde_matches = mol.GetSubstructMatches(aromatic_aldehyde_pattern)
    n_aromatic_aldehydes = len(aromatic_aldehyde_matches)
    
    if n_aromatic_aldehydes > 0:
        return True, f"Contains {n_aromatic_aldehydes} aldehyde group(s) attached to an aromatic moiety out of {n_aldehydes} aldehyde group(s) in the molecule"
    else:
        return False, f"No aldehyde group attached to an aromatic moiety found among {n_aldehydes} aldehyde group(s) in the molecule"