"""
Classifies: CHEBI:26377 pterocarpans
"""
from rdkit import Chem

def is_pterocarpans(smiles: str):
    """
    Determines if a molecule is a part of the pterocarpan class based on its SMILES string.
    A pterocarpan is characterized by the 6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a pterocarpan, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern to represent the pterocarpan skeleton
    # This pattern focuses on the tricyclic core with the right connectivity for a pterocarpan skeleton
    pterocarpan_pattern = Chem.MolFromSmarts("C1Oc2ccccc2C3Cc4ccccc4OC13")

    # Check if the molecule has the pterocarpan skeleton
    if mol.HasSubstructMatch(pterocarpan_pattern):
        return True, "Contains 6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene skeleton"
    else:
        return False, "Does not contain the defining pterocarpan skeleton"