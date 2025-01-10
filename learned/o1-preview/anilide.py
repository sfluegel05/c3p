"""
Classifies: CHEBI:13248 anilide
"""
from rdkit import Chem

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.
    An anilide is an aromatic amide obtained by acylation of aniline,
    meaning it contains an amide group where the nitrogen atom is directly attached to an aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anilide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the anilide SMARTS pattern:
    # Aromatic carbon ([c]) connected to nitrogen ([N]) connected to carbonyl carbon (C=O)
    anilide_pattern = Chem.MolFromSmarts("[c][NH0,H1][C](=O)[#6]")

    if mol.HasSubstructMatch(anilide_pattern):
        return True, "Contains an anilide functional group (amide with nitrogen attached to an aromatic ring)"
    else:
        return False, "Anilide substructure not found"