"""
Classifies: CHEBI:39437 tocol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tocol(smiles: str):
    """
    Determines if a molecule is a tocol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a tocol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Broadly define chroman/chromene core structure: aromatic + oxygen in six-membered ring
    chroman_core_pattern = Chem.MolFromSmarts("C1=CC(O)C=C1")

    if not mol.HasSubstructMatch(chroman_core_pattern):
        return False, "No chroman or chromene moiety detected"

    # Check for a sufficient hydrocarbon tail that typically characterizes tocols.
    long_hydrocarbon_pattern = Chem.MolFromSmarts("C(C)C")
    tail_matches = mol.GetSubstructMatches(long_hydrocarbon_pattern)
    if len(tail_matches) < 3:
        return False, "Insufficient hydrocarbon tail length"

    # Make sure at least one hydroxyl group is attached to the aromatic core
    hydroxyl_pattern = Chem.MolFromSmarts("cO")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No aromatic-linked hydroxyl found"

    return True, "Contains appropriate chroman/chromene moiety with hydrocarbon tail and hydroxyl group"