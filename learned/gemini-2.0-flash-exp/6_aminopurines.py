"""
Classifies: CHEBI:20706 6-aminopurines
"""
from rdkit import Chem

def is_6_aminopurines(smiles: str):
    """
    Determines if a molecule contains a 6-aminopurine (adenine) moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains 6-aminopurine, False otherwise
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for the 6-aminopurine (adenine) moiety
    adenine_smarts = Chem.MolFromSmarts("n1cnc2c(n1)nc([NH2])n2")

    if adenine_smarts is None:
        return False, "Invalid SMARTS string"

    # Check for the presence of the 6-aminopurine substructure
    if mol.HasSubstructMatch(adenine_smarts):
        return True, "Contains a 6-aminopurine moiety"
    else:
        return False, "Does not contain a 6-aminopurine moiety"