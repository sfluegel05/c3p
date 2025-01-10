"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
from rdkit import Chem

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.
    An aliphatic alcohol is defined as an alcohol derived from an aliphatic compound,
    with at least one hydroxyl group attached to an aliphatic carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aliphatic alcohol, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts('[CX4][OH]')
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group attached to an aliphatic carbon found"

    # Check if all carbons are non-aromatic (aliphatic)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            if atom.GetIsAromatic():
                return False, "Contains aromatic carbon atoms"

    return True, "Molecule contains a hydroxyl group attached to an aliphatic carbon, indicative of an aliphatic alcohol"