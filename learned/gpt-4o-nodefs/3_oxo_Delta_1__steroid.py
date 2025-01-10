"""
Classifies: CHEBI:20156 3-oxo-Delta(1) steroid
"""
from rdkit import Chem

def is_3_oxo_Delta_1_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(1) steroid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-Delta(1) steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a broad steroid backbone pattern with a 3-keto group and Delta(1) double bond
    # The carbonyl group (C=O) at the third position and a double bond between the first and second carbon.
    steroid_pattern = Chem.MolFromSmarts("C1=CC2CC[C@H]3C(=O)CC[C@@]3(C)C2CC1")

    # Check for the steroid pattern in the molecule
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "3-oxo-Delta(1) steroid structure not found"

    return True, "Molecule contains a 3-oxo-Delta(1) steroid structure"