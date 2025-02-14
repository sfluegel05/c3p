"""
Classifies: CHEBI:140310 phenyl acetates
"""
"""
Classifies: Phenyl Acetates
"""
from rdkit import Chem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate based on its SMILES string.
    A phenyl acetate is defined as an acetate ester obtained by formal condensation 
    of the carboxy group of acetic acid with the hydroxy group of any phenol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a phenyl acetate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for phenyl acetate group
    phenyl_acetate_smarts = '[CH3][C](=O)[O][c]'
    phenyl_acetate_pattern = Chem.MolFromSmarts(phenyl_acetate_smarts)
    if phenyl_acetate_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Search for phenyl acetate pattern in the molecule
    matches = mol.GetSubstructMatches(phenyl_acetate_pattern)
    if matches:
        return True, "Contains phenyl acetate group"
    else:
        return False, "No phenyl acetate group found"