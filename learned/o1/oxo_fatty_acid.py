"""
Classifies: CHEBI:59644 oxo fatty acid
"""
from rdkit import Chem

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.
    An oxo fatty acid is defined as any fatty acid containing at least one aldehydic or ketonic group
    in addition to the carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oxo fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify carboxylic acid groups (-C(=O)OH)
    carboxylic_acid = Chem.MolFromSmarts("C(=O)[O;H1]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid)
    if not carboxylic_acid_matches:
        return False, "No carboxylic acid group found"

    # Identify aldehyde groups (-CHO)
    aldehyde = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde)

    # Identify ketone groups (>C=O)
    ketone = Chem.MolFromSmarts("[#6][CX3](=O)[#6]")
    ketone_matches = mol.GetSubstructMatches(ketone)

    # Exclude the carboxylic acid group from matches
    carboxylic_acid_atoms = set()
    for match in carboxylic_acid_matches:
        carboxylic_acid_atoms.update(match)

    # Filter out oxo groups that are part of the carboxylic acid
    oxo_matches = []
    for match in aldehyde_matches + ketone_matches:
        if not carboxylic_acid_atoms.intersection(match):
            oxo_matches.append(match)

    if not oxo_matches:
        return False, "No additional aldehydic or ketonic group found"

    # Check for sufficient chain length to be a fatty acid (minimum 4 carbons)
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbon_atoms) < 4:
        return False, f"Insufficient carbon atoms ({len(carbon_atoms)}) for a fatty acid"

    return True, "Contains carboxylic acid and additional oxo group; classified as an oxo fatty acid"