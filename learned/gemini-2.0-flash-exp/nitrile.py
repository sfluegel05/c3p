"""
Classifies: CHEBI:18379 nitrile
"""
from rdkit import Chem

def is_nitrile(smiles: str):
    """
    Determines if a molecule is a nitrile based on its SMILES string.
    A nitrile is characterized by a carbon atom triply bonded to a nitrogen atom (C#N).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrile, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a nitrile group (C#N)
    nitrile_pattern = Chem.MolFromSmarts("C#N")

    # Check if the molecule has the nitrile substructure
    if mol.HasSubstructMatch(nitrile_pattern):
        return True, "Contains a nitrile group (C#N)"
    else:
        return False, "Does not contain a nitrile group (C#N)"