"""
Classifies: CHEBI:134363 tertiary amine oxide
"""
from rdkit import Chem

def is_tertiary_amine_oxide(smiles: str):
    """
    Determines if a molecule is a tertiary amine oxide based on its SMILES string.
    A tertiary amine oxide consists of an N-oxide where there are three organic groups bonded to the nitrogen atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine oxide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Extended pattern to match an N+ bonded to any three organic groups (considering a variety of substitutions),
    # and a single negatively charged oxygen (N-O) ensuring flexibility.
    tertiary_amine_oxide_pattern = Chem.MolFromSmarts("[NX4+]([#6,#1,#7,#8,#9,#17,#35,#53])([#6,#1,#7,#8,#9,#17,#35,#53])([#6,#1,#7,#8,#9,#17,#35,#53])[O-]")

    if mol.HasSubstructMatch(tertiary_amine_oxide_pattern):
        return True, "Molecule contains a tertiary amine oxide group"

    return False, "Molecule does not contain a tertiary amine oxide group"