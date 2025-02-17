"""
Classifies: CHEBI:52221 isothiocyanate
"""
"""
Classifies: Isothiocyanate (R-N=C=S, an organosulfur compound)
"""
from rdkit import Chem

def is_isothiocyanate(smiles: str):
    """
    Determines if a molecule belongs to the isothiocyanate class based on its SMILES string.
    An isothiocyanate is defined as an organosulfur compound with the general formula R-N=C=S.
    
    This function checks for the presence of the characteristic isothiocyanate group.
    Note: Some molecules may present the isothiocyanate group as S=C=N. Both representations are considered.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an isothiocyanate, False otherwise.
        str: Reason for the classification.
    """
    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for isothiocyanate.
    # Pattern 1: as given in the definition R-N=C=S
    pattern1 = Chem.MolFromSmarts("N=C=S")
    # Pattern 2: sometimes the group is represented as S=C=N in the SMILES
    pattern2 = Chem.MolFromSmarts("S=C=N")
    
    # Check if the molecular structure contains either pattern
    if mol.HasSubstructMatch(pattern1) or mol.HasSubstructMatch(pattern2):
        return True, "Contains isothiocyanate group (R-N=C=S)"
    else:
        return False, "Does not contain an isothiocyanate group"