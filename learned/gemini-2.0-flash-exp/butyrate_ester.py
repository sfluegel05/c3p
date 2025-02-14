"""
Classifies: CHEBI:50477 butyrate ester
"""
from rdkit import Chem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.
    A butyrate ester contains a butyryl group (CH3CH2CH2C(=O)-) linked to an oxygen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butyrate ester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the butyrate ester substructure using SMARTS pattern.
    butyrate_pattern = Chem.MolFromSmarts("CCCC(=O)O")

    if mol.HasSubstructMatch(butyrate_pattern):
        return True, "Molecule contains a butyrate ester group"
    else:
        return False, "No butyrate ester substructure found"