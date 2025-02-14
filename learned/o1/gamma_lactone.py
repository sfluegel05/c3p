"""
Classifies: CHEBI:37581 gamma-lactone
"""
"""
Classifies: gamma-lactone
"""

from rdkit import Chem

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule is a gamma-lactone based on its SMILES string.
    A gamma-lactone is defined as a lactone (cyclic ester) having a five-membered ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a gamma-lactone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string to RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for gamma-lactone
    gamma_lactone_smarts = '[#6;R]=O-[O;R]-[#6;R]-[#6;R]-[#6;R]'
    gamma_lactone_pattern = Chem.MolFromSmarts(gamma_lactone_smarts)

    # Search for gamma-lactone pattern in the molecule
    matches = mol.GetSubstructMatches(gamma_lactone_pattern)
    if matches:
        return True, "Contains gamma-lactone ring (five-membered lactone ring)"
    else:
        return False, "No gamma-lactone ring found"