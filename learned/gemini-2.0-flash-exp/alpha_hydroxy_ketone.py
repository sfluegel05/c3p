"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.
    An alpha-hydroxy ketone is a ketone with a hydroxy group on the alpha-carbon relative to the C=O group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    # 1. Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 2. Define SMARTS pattern for alpha-hydroxy ketone
    alpha_hydroxy_ketone_pattern = Chem.MolFromSmarts("[CX3](=O)[CX4H1,CX4H2][OX2H1]") 
    
    # 3. Check for the pattern
    matches = mol.GetSubstructMatches(alpha_hydroxy_ketone_pattern)

    if matches:
      return True, "Hydroxy group found on an alpha carbon of a ketone"
    else:
        return False, "No hydroxy group found on any alpha carbon of a ketone"