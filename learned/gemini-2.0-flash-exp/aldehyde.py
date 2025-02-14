"""
Classifies: CHEBI:17478 aldehyde
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aldehyde(smiles: str):
    """
    Determines if a molecule is an aldehyde based on its SMILES string.
    An aldehyde is a compound with a carbonyl group (C=O) bonded to a hydrogen (or deuterium) and another R group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldehyde, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Create SMARTS pattern for aldehyde group (C=O-H or C=O-D, and one R group).
    aldehyde_pattern = Chem.MolFromSmarts("[CX3](=O)[H,D][!H]")

    # Check for the presence of the aldehyde group
    if mol.HasSubstructMatch(aldehyde_pattern):
        return True, "Molecule contains an aldehyde group (C=O-H or C=O-D, and one R group)"
    else:
        return False, "Molecule does not contain an aldehyde group (C=O-H or C=O-D, and one R group)"