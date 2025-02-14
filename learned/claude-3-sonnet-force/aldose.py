"""
Classifies: CHEBI:15693 aldose
"""
"""
Classifies: CHEBI:17395 aldose
An aldose is an aldehydic parent sugar (polyhydroxy aldehyde H[CH(OH)]nC(=O)H, n >= 2)
or its intramolecular hemiacetal.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_aldose(smiles: str):
    """
    Determines if a molecule is an aldose based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldose, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for aldehyde or hemiacetal group
    aldehyde_pattern = Chem.MolFromSmarts("[CH](=O)[OH]")
    hemiacetal_pattern = Chem.MolFromSmarts("[CH1](O)[OH]")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    hemiacetal_matches = mol.GetSubstructMatches(hemiacetal_pattern)
    if not aldehyde_matches and not hemiacetal_matches:
        return False, "No aldehyde or hemiacetal group found"

    # Check for multiple hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    n_hydroxyls = len(hydroxyl_matches)
    if n_hydroxyls < 2:
        return False, "Less than 2 hydroxyl groups found"

    # Count carbon atoms
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if n_carbons < 3:
        return False, "Less than 3 carbon atoms found"

    # Check for incompatible functional groups
    incompatible_patterns = [
        Chem.MolFromSmarts("C(=O)O"),  # Carboxylic acid
        Chem.MolFromSmarts("N"),  # Amine
        Chem.MolFromSmarts("O=C-O-C"),  # Ester
        Chem.MolFromSmarts("C-O-C"),  # Ether
    ]
    for pattern in incompatible_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, f"Incompatible functional group found: {pattern.GetSmarts()}"

    # Passed all checks
    return True, "Molecule meets structural requirements for an aldose"