"""
Classifies: CHEBI:35757 monocarboxylic acid anion
"""
"""
Classifies: CHEBI:24868 monocarboxylic acid anion
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a monocarboxylic acid anion based on its SMILES string.
    A monocarboxylic acid anion is formed when the carboxy group of a monocarboxylic acid is deprotonated.

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

    # Look for carboxylate anion pattern (-C(=O)[O-])
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)

    # Check if there is exactly one carboxylate group
    if len(carboxylate_matches) != 1:
        return False, f"Found {len(carboxylate_matches)} carboxylate groups, should be exactly 1"

    # Count number of carboxy groups (-C(=O)OH)
    carboxy_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxy_matches = mol.GetSubstructMatches(carboxy_pattern)

    # Ensure there are no carboxy groups (i.e., all have been deprotonated)
    if len(carboxy_matches) > 0:
        return False, "Found intact carboxy groups, should be deprotonated"

    # Check for multiple carboxylate groups
    dicarboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-].[C@@](C(=O)[O-])(C)C")
    if mol.HasSubstructMatch(dicarboxylate_pattern):
        return False, "Contains multiple carboxylate groups"

    # Count number of ring systems
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()

    # Monocarboxylic acid anions typically contain at most one ring system
    if num_rings > 1:
        return False, "Contains more than one ring system"

    return True, "Contains a single carboxylate anion group"