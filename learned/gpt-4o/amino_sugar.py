"""
Classifies: CHEBI:28963 amino sugar
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    An amino sugar is characterized by one or more amino groups replacing hydroxyl groups on a sugar molecule.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino sugar, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS for amino group which should replace a hydroxyl group
    amino_group_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")  # Primary or secondary amine not in amides

    # Check we have an amino group
    if not mol.HasSubstructMatch(amino_group_pattern):
        return False, "No amino groups found"

    # Basic pattern for a monosaccharide ring (e.g., hexopyranose)
    sugar_ring_pattern = Chem.MolFromSmarts("C1(OC)[C@H](O)[C@@H](O)[C@H](O)[C@H](O)C[O]1")

    # Check for sugar-like structure
    if not mol.HasSubstructMatch(sugar_ring_pattern):
        return False, "No sugar-like structure found"

    # Analyze if the amino group replaces an OH: We can check this by ensuring no free hydroxyl on same carbon
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N':
            # Inspect bonded atoms to check if one of the oxygens (replaced hydroxyl) is missing
            neighbors = [a.GetSymbol() for a in atom.GetNeighbors()]
            if 'C' in neighbors:
                carbon = next(a for a in atom.GetNeighbors() if a.GetSymbol() == 'C')
                _oh_neighbors = [a.GetSymbol() for a in carbon.GetNeighbors()]
                if 'O' in _oh_neighbors:
                    continue  # Still has oxygen, so not replacing the hydroxyl group
                return True, "Contains sugar backbone with one or more amino groups replacing hydroxy groups"

    return False, "No correct hydroxyl replacement was found by an amino group"