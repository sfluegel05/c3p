"""
Classifies: CHEBI:134251 guaiacols
"""
from rdkit import Chem

def is_guaiacols(smiles: str):
    """
    Determines if a molecule is a guaiacol based on its SMILES string.
    Guaiacols are characterized by a specific arrangement of a methoxy group and a hydroxyl group on a benzene ring.
    
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

    # Guaiacol core pattern: A benzene ring with a methoxy (OC) group ortho to a hydroxyl (O) group.
    # Example search pattern that is more specific than the general aromatic where additional hydrogens are considered.
    guaiacol_pattern = Chem.MolFromSmarts("c1cc(OC)c(O)cc1")  # This pattern previously caused false positives
    
    # Introducing a refined pattern, ensuring the structure represents typical guaiacol more strictly
    refined_guaiacol_pattern = Chem.MolFromSmarts("c1c(O)cc(OC)cc1")  # Refined pattern with correct ortho orientation
    
    # Apply refined guaiacol SMARTS pattern matching
    if mol.HasSubstructMatch(refined_guaiacol_pattern):
        return True, "Contains guaiacol core structure with refined pattern"

    return False, "No guaiacol core structure found with refined pattern"