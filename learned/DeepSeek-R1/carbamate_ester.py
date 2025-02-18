"""
Classifies: CHEBI:23003 carbamate ester
"""
"""
Classifies: CHEBI:23004 carbamate ester
"""
from rdkit import Chem

def is_carbamate_ester(smiles: str):
    """
    Determines if a molecule is a carbamate ester based on its SMILES string.
    A carbamate ester is an ester of carbamic acid (structure: R-O-C(=O)-N<).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbamate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define carbamate ester pattern: [C]-O-C(=O)-N where N has valence 3
    # Ensures oxygen is part of ester group (connected to R and C=O)
    # [N;v3] matches nitrogens with valence 3 (including aromatic)
    carbamate_pattern = Chem.MolFromSmarts("[C][OX2]C(=O)[N;v3]")

    # Check for the presence of the carbamate group
    if mol.HasSubstructMatch(carbamate_pattern):
        return True, "Contains carbamate ester group (O-C(=O)-N)"
    else:
        return False, "No carbamate ester group found"