"""
Classifies: CHEBI:35757 monocarboxylic acid anion
"""
"""
Classifies: CHEBI:35537 monocarboxylic acid anion

A carboxylic acid anion formed when the carboxy group of a monocarboxylic acid is deprotonated.
"""

from rdkit import Chem

def is_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a monocarboxylic acid anion based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monocarboxylic acid anion, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of carboxylate anion (-COO-)
    carboxylate_pattern = Chem.MolFromSmarts("[O-]C(=O)")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate anion (-COO-) found"

    # Check for exactly one carboxylate anion
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    if len(carboxylate_matches) != 1:
        return False, f"Found {len(carboxylate_matches)} carboxylate anions, expected exactly 1"

    # Check that there is only one carboxylate group
    carbonyl_oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == 0)
    if carbonyl_oxygen_count > 1:
        return False, "Found multiple carbonyl oxygens, expected a monocarboxylic acid"

    return True, "Contains a single carboxylate anion (-COO-) group"