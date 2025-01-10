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

    # Refined pattern to match a nitrogen with a positive charge bonded to three carbons and one oxygen with negatively charged oxygen (N+ with C and O-)
    tertiary_amine_oxide_pattern = Chem.MolFromSmarts("[N+X4]([CX4,a])([CX4,a])([CX4,a])[OX1-]")
    
    if mol.HasSubstructMatch(tertiary_amine_oxide_pattern):
        return True, "Molecule contains a tertiary amine oxide group"

    return False, "Molecule does not contain a tertiary amine oxide group"