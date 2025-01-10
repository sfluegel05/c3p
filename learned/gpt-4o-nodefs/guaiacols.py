"""
Classifies: CHEBI:134251 guaiacols
"""
from rdkit import Chem

def is_guaiacols(smiles: str):
    """
    Determines if a molecule is a guaiacol based on its SMILES string.
    Guaiacols are characterized by a methoxy group and a hydroxyl group on a benzene ring.
    This function tries to detect these groups' presence in characteristic positions uniquely.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a guaiacol, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Guaiacol core patterns: benzene ring with hydroxyl and methoxy adjacent positions.
    # We ensure positions usually found in guaiacols, for example, position ('O') next to benzene rings.
    guaiacol_pattern = Chem.MolFromSmarts("c1cc(OC)c(O)cc1")
    
    if mol.HasSubstructMatch(guaiacol_pattern):
        return True, "Contains guaiacol core structure"
    
    return False, "No guaiacol core structure found"

# Additional patterns can be added here if more diverse structures of guaiacols are to be included.