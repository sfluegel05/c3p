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

    # Check for the carboxylic acid group (-C(=O)OH)
    carboxylic_acid = Chem.MolFromSmarts("C(=O)[O;H1]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid)
    if not carboxylic_acid_matches:
        return False, "No carboxylic acid group found"

    # Get the atoms involved in carboxylic acid groups
    carboxylic_acid_atoms = set()
    for match in carboxylic_acid_matches:
        carboxylic_acid_atoms.update(match)

    # Check for aldehyde groups (-CHO)
    aldehyde = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    aldehyde_matches = []
    for match in mol.GetSubstructMatches(aldehyde):
        if set(match).isdisjoint(carboxylic_acid_atoms):
            aldehyde_matches.append(match)

    # Check for ketone groups (>C=O)
    ketone = Chem.MolFromSmarts("[#6][CX3](=O)[#6]")
    ketone_matches = []
    for match in mol.GetSubstructMatches(ketone):
        if set(match).isdisjoint(carboxylic_acid_atoms):
            ketone_matches.append(match)

    total_oxo_groups = len(aldehyde_matches) + len(ketone_matches)

    if total_oxo_groups == 0:
        return False, "No additional aldehydic or ketonic group found"

    # Check for aliphatic chain length (minimum 4 carbons)
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(c_atoms) < 4:
        return False, "Aliphatic chain is too short to be a fatty acid"

    return True, "Contains carboxylic acid and additional aldehydic or ketonic group"