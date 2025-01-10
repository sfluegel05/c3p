"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
from rdkit import Chem

def is_7_hydroxyisoflavones(smiles: str):
    """
    Determines if a molecule is a 7-hydroxyisoflavone based on its SMILES string.
    A 7-hydroxyisoflavone possesses a hydroxyisoflavone skeleton with a hydroxy group specifically at the 7-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 7-hydroxyisoflavone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define core isoflavone skeleton (benzopyran core with attached phenyl ring and carbonyl)
    # This pattern needs to be flexible to match various substitutions and ring closures
    # Key feature: pyran ring with adjacent aromatic and carbonyl group
    isoflavone_core_pattern = Chem.MolFromSmarts("O=C1C2=COC3=CC=CC=C3C2=CC=C1")

    if not mol.HasSubstructMatch(isoflavone_core_pattern):
        return False, "No isoflavone core structure found"

    # Flexible pattern that identifies hydroxy group at the 7-position in the simplest A-ring scenario
    hydroxy_pattern = Chem.MolFromSmarts("O-c1cc2occc2c(=O)c1")  # Hydroxy group at C7 position of A ring

    # Check for 7-hydroxy position specifically
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No hydroxy group at the 7-position found"

    return True, "Has isoflavone skeleton with hydroxy group at 7-position"