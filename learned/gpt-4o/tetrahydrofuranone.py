"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
from rdkit import Chem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    A tetrahydrofuranone is any oxolane with an oxo-substituent on the tetrahydrofuran ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrahydrofuranone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Define a SMARTS pattern for a five-membered ring with an oxygen and an oxo group at any position
    thf_with_oxo = Chem.MolFromSmarts("O1[C;R][C;R][C;R](=O)[C;R]1 |or1:O1[C;R](=O)[C;R][C;R][C;R]1|")

    if thf_with_oxo is None:
        return None, "Error in generating SMARTS pattern"

    # Check against the pattern
    if mol.HasSubstructMatch(thf_with_oxo):
        return True, "Contains a tetrahydrofuran ring with an oxo group"

    return False, "Does not contain a tetrahydrofuran ring with an oxo group"