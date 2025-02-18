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

    # Define carbamate ester pattern: [O connected to C]-C(=O)-N
    # The oxygen must be part of an ester group (connected to two carbons)
    carbamate_pattern = Chem.MolFromSmarts("[OX2]C(=O)[NX3]")

    # Check for the presence of the carbamate group
    if mol.HasSubstructMatch(carbamate_pattern):
        return True, "Contains carbamate ester group (O-C(=O)-N)"
    else:
        return False, "No carbamate ester group found"