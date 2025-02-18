"""
Classifies: CHEBI:22307 aldoxime
"""
"""
Classifies: CHEBI:36576 aldoxime
"""
from rdkit import Chem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime based on its SMILES string.
    Aldoximes are oximes of aldehydes with the general structure RCH=NOH.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldoxime, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the aldoxime SMARTS pattern: [C;D1,D2]=N-O
    # Matches carbons with degree 1 or 2 (from aldehyde) connected to oxime group
    aldoxime_pattern = Chem.MolFromSmarts('[C;D1,D2]=[NX2]-[OX2H]')
    
    # Check for substructure match
    if mol.HasSubstructMatch(aldoxime_pattern):
        return True, "Contains aldoxime group (R-CH=N-OH)"
    else:
        return False, "No aldoxime group detected"