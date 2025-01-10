"""
Classifies: CHEBI:134251 guaiacols
"""
"""
Classifies: Guaiacols
"""
from rdkit import Chem

def is_guaiacols(smiles: str):
    """
    Determines if a molecule is a guaiacol based on its SMILES string.
    Guaiacols typically have a methoxy group and a hydroxyl group attached to a benzene ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a guaiacol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the guaiacol pattern: a benzene ring with adjacent methoxy and hydroxyl groups
    guaiacol_pattern = Chem.MolFromSmarts("c1cc(OC)c(O)cc1")
    if mol.HasSubstructMatch(guaiacol_pattern):
        return True, "Contains guaiacol core structure"

    return False, "No guaiacol core structure found"