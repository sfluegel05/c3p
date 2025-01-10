"""
Classifies: CHEBI:22307 aldoxime
"""
from rdkit import Chem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime based on its SMILES string.
    An aldoxime contains a carbon-nitrogen double bond (C=N) and a hydroxyl group
    directly attached to the nitrogen (C=N-O structure), where the carbon is part
    of an aldehyde group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aldoxime, False otherwise
        str: Reason for classification
    """

    # Parse SMILES into RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for an aldoxime, ensuring the carbon is part of an aldehyde
    # Carbon with one single bond (e.g., C=NO with hydrogen attachment)
    # However, we cannot strictly enforce only one other attachment due to RS conformation variability.
    aldoxime_pattern = Chem.MolFromSmarts("[CX3]=[NX2]O")  # '[CX3]=[NX2]O' better reflects stricter structure

    if mol.HasSubstructMatch(aldoxime_pattern):
        return True, "Contains the characteristic C=N-O structure of aldoximes"

    return False, "Does not contain the characteristic C=N-O structure of aldoximes"