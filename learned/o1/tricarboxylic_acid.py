"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
"""
Classifies: tricarboxylic acid
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid based on its SMILES string.
    A tricarboxylic acid is an oxoacid containing exactly three carboxy groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tricarboxylic acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Sanitize molecule
    Chem.SanitizeMol(mol)

    # Identify carboxy groups (both protonated and deprotonated forms)
    carboxy_pattern = Chem.MolFromSmarts('[CX3](=O)[O-,OH]')
    carboxy_matches = mol.GetSubstructMatches(carboxy_pattern)
    num_carboxy_groups = len(carboxy_matches)

    if num_carboxy_groups != 3:
        return False, f"Found {num_carboxy_groups} carboxy groups, need exactly 3"

    # Check for other acidic groups (e.g., sulfonic acids, phosphonic acids)
    other_acid_patterns = [
        Chem.MolFromSmarts('[SX4](=O)(=O)[O-,OH]'),  # sulfonic acid
        Chem.MolFromSmarts('[PX4](=O)([O-,OH])[O-,OH]'),  # phosphonic acid
    ]
    for pattern in other_acid_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, "Contains other acidic groups (e.g., sulfonic or phosphonic acids)"

    # Ensure molecule is organic (contains carbon)
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Molecule does not contain carbon, not an organic acid"

    return True, "Contains exactly three carboxy groups (tricarboxylic acid)"