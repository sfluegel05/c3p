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

    # Define the porphyrin core as a SMARTS pattern
    porphyrin_smarts = 'c1cc2cc3ccc(cc4ccc(cc1n2)n4)n3'  # Porphyrin core without hydrogen specification
    porphyrin_pattern = Chem.MolFromSmarts(porphyrin_smarts)

    if porphyrin_pattern is None:
        return False, "Invalid porphyrin SMARTS pattern"

    # Perform substructure search for the porphyrin core
    if mol.HasSubstructMatch(porphyrin_pattern):
        return True, "Porphyrin ring system detected"
    else:
        return False, "No porphyrin ring system detected"