"""
Classifies: CHEBI:50477 butyrate ester
"""
from rdkit import Chem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.
    A butyrate ester contains a butyric acid component esterified with an alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butyrate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more comprehensive pattern for butyrate ester
    # Target 'CCC(=O)O' for the butyric acid part with flexibly integrated esters (alkyl/aryl)
    # Extend search to match potential substituent variations and longer tail chains
    butyrate_pattern = Chem.MolFromSmarts("[$([CH3]C(=O)O),$([CH2]CC(=O)O),$([CH2][CH2]C(=O)O)]")

    # Check for presence of butyric acid ester component
    if mol.HasSubstructMatch(butyrate_pattern):
        return True, "Contains butyric acid ester linkage"

    return False, "No butyric acid ester linkage found"