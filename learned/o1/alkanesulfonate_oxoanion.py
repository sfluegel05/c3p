"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
"""
Classifies: alkanesulfonate oxoanion
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    An alkanesulfonate oxoanion is characterized by a sulfonate group (-S(=O)(=O)[O-])
    attached to a carbon atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanesulfonate oxoanion, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the alkanesulfonate oxoanion pattern
    alkanesulfonate_pattern = Chem.MolFromSmarts("[#6]-[S](=O)(=O)[O-]")
    if alkanesulfonate_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Search for the pattern in the molecule
    matches = mol.GetSubstructMatches(alkanesulfonate_pattern)
    if not matches:
        return False, "Does not contain alkanesulfonate oxoanion group"

    # Verify that the sulfur is connected to a carbon atom
    for match in matches:
        sulfur_idx = match[1]  # Index of the sulfur atom in the match
        sulfur_atom = mol.GetAtomWithIdx(sulfur_idx)
        neighbors = [atom.GetAtomicNum() for atom in sulfur_atom.GetNeighbors()]
        if 6 in neighbors:
            return True, "Contains alkanesulfonate oxoanion group attached to carbon"

    return False, "Sulfonate group not attached to carbon atom"