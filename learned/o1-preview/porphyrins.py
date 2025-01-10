"""
Classifies: CHEBI:26214 porphyrins
"""
"""
Classifies: CHEBI:8338 porphyrins
"""

from rdkit import Chem

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    A porphyrin contains a fundamental skeleton of four pyrrole nuclei united through the alpha-positions by four methine groups to form a macrocyclic structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a porphyrin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the porphyrin core using SMARTS pattern
    porphyrin_smarts = """
    [$([nH]),N]1C=CC2=NC=CC3=NC=CC4=NC=CC=1C2=CC3=CC4
    """

    # Remove newlines and whitespace from SMARTS
    porphyrin_smarts = porphyrin_smarts.replace('\n', '').replace(' ', '')
    porphyrin_pattern = Chem.MolFromSmarts(porphyrin_smarts)
    if porphyrin_pattern is None:
        return False, "Error in porphyrin SMARTS pattern"

    # Perform substructure search for porphyrin core
    if mol.HasSubstructMatch(porphyrin_pattern):
        return True, "Porphyrin core structure detected"
    else:
        return False, "No porphyrin core structure detected"