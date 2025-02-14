"""
Classifies: CHEBI:26214 porphyrins
"""
"""
Classifies: Porphyrins
"""
from rdkit import Chem

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    A porphyrin is characterized by a macrocyclic structure composed of four pyrrole rings connected via methine bridges.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a porphyrin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the porphyrin core SMARTS pattern
    # The pattern represents the macrocyclic structure of porphyrins
    porphyrin_smarts = 'n1c(ccc1)-c1c2c(ccc2)[nH]c2ccc(-c3ccc([nH]3)c(-c3ccc([nH]3)c1))c2'
    porphyrin_pattern = Chem.MolFromSmarts(porphyrin_smarts)
    if porphyrin_pattern is None:
        return None, "Failed to create porphyrin SMARTS pattern"

    # Check for porphyrin core
    if mol.HasSubstructMatch(porphyrin_pattern):
        return True, "Contains porphyrin core structure with four pyrrole rings connected via methine bridges"
    else:
        return False, "Does not contain the porphyrin core macrocyclic structure"