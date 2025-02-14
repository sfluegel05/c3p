"""
Classifies: CHEBI:26377 pterocarpans
"""
from rdkit import Chem

def is_pterocarpans(smiles: str):
    """
    Determines if a molecule is part of the pterocarpan class based on its SMILES string.
    A pterocarpan has a 6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pterocarpan, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more general pterocarpan skeleton SMARTS pattern
    pterocarpan_pattern = Chem.MolFromSmarts("C12Oc3ccccc3C(C4C=CC=CC4)=C1COC2")

    # Check for the presence of the pterocarpan skeleton
    if mol.HasSubstructMatch(pterocarpan_pattern):
        return True, "Contains 6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene skeleton"
    else:
        return False, "Does not contain the defining pterocarpan skeleton"