"""
Classifies: CHEBI:65061 2,5-diketopiperazines
"""
"""
Classifies: CHEBI:27206 2,5-diketopiperazine
"""
from rdkit import Chem

def is_2_5_diketopiperazines(smiles: str):
    """
    Determines if a molecule is a 2,5-diketopiperazine based on its SMILES string.
    A 2,5-diketopiperazine is a six-membered ring containing two amide (-C(=O)N-) linkages at positions 1 and 4, with carbonyl groups at positions 2 and 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2,5-diketopiperazine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a 2,5-diketopiperazine core. Explicitly specifies a 6 membered ring.
    diketopiperazine_pattern = Chem.MolFromSmarts("[N1]C(=O)[C][N2]C(=O)[C]12")
    
    # Check if the molecule matches the 2,5-diketopiperazine pattern
    if mol.HasSubstructMatch(diketopiperazine_pattern):
       return True, "Molecule contains a 2,5-diketopiperazine core structure"
    else:
        return False, "Molecule does not contain a 2,5-diketopiperazine core structure"