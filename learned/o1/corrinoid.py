"""
Classifies: CHEBI:33913 corrinoid
"""
"""
Classifies: CHEBI:33909 corrinoid
"""

from rdkit import Chem

def is_corrinoid(smiles: str):
    """
    Determines if a molecule is a corrinoid based on its SMILES string.
    A corrinoid is a derivative of the corrin nucleus, which contains four reduced or partly reduced
    pyrrole rings joined in a macrocycle by three =C- groups and one direct carbon-carbon bond
    linking alpha positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a corrinoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the corrin nucleus as a SMARTS pattern
    corrin_smarts = '[nH]1ccc2c1ccc1c(c2)ccc2c1[nH]ccc2'  # Simplified corrin core

    pattern = Chem.MolFromSmarts(corrin_smarts)
    if pattern is None:
        return False, "Invalid corrin nucleus SMARTS pattern"

    # Check if the molecule contains the corrin nucleus
    if mol.HasSubstructMatch(pattern):
        return True, "Molecule contains the corrin nucleus"
    else:
        return False, "Molecule does not contain the corrin nucleus"