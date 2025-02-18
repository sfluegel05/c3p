"""
Classifies: CHEBI:50753 isoflavonoid
"""
"""
Classifies: CHEBI:18214 isoflavonoid
"""
from rdkit import Chem

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    An isoflavonoid is a 1-benzopyran (chromene) with an aryl substituent at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the core chromene (1-benzopyran) structure with ketone at position 4
    chromene_core = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[O][C](=[O])[C]=[C]")
    if not mol.HasSubstructMatch(chromene_core):
        return False, "No chromene (1-benzopyran) core found"
    
    # Check for aryl substituent at position 3 (carbon adjacent to oxygen in the pyran ring)
    aryl_pattern = Chem.MolFromSmarts("[O][C](=[O])[C]([c;R])=[C]")
    if not mol.HasSubstructMatch(aryl_pattern):
        return False, "No aryl substituent at position 3"
    
    return True, "Contains chromene core with aryl group at position 3"