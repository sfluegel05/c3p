"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
"""
Classifies: omega-hydroxy fatty acid
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an omega-hydroxy fatty acid based on its SMILES string.
    An omega-hydroxy fatty acid is a straight-chain fatty acid with a carboxyl group
    at position 1 and a hydroxyl group at the omega position (opposite end).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an omega-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify carboxyl group [C(=O)O]
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1H]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No carboxyl group found"

    # Identify hydroxyl groups [O-H]
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if not hydroxyl_matches:
        return False, "No hydroxyl group found"

    # Get the indices of the carboxyl carbon(s) and hydroxyl oxygen(s)
    carboxyl_c_indices = [match[0] for match in carboxyl_matches]
    hydroxyl_o_indices = [match[0] for match in hydroxyl_matches]

    # Check if carboxyl group is at terminal position
    is_terminal_carboxyl = False
    for c_idx in carboxyl_c_indices:
        atom = mol.GetAtomWithIdx(c_idx)
        # Carboxyl carbon should be connected to one heavy atom (excluding O)
        neighbors = [a for a in atom.GetNeighbors() if a.GetAtomicNum() > 1 and a.GetAtomicNum() != 8]
        if len(neighbors) == 1:
            carboxyl_c_idx = c_idx
            is_terminal_carboxyl = True
            break
    if not is_terminal_carboxyl:
        return False, "Carboxyl group is not at terminal position"

    # Check if hydroxyl group is at terminal position
    is_terminal_hydroxyl = False
    for o_idx in hydroxyl_o_indices:
        oxygen = mol.GetAtomWithIdx(o_idx)
        # Hydroxyl oxygen should be connected to one heavy atom (the carbon)
        if len(oxygen.GetNeighbors()) != 1:
            continue
        carbon = oxygen.GetNeighbors()[0]
        if carbon.GetAtomicNum() == 6:
            # Check if carbon is terminal (connected to only one heavy atom)
            heavy_atom_neighbors = [a for a in carbon.GetNeighbors() if a.GetAtomicNum() > 1]
            if len(heavy_atom_neighbors) == 1:
                hydroxyl_o_idx = o_idx
                hydroxyl_c_idx = carbon.GetIdx()
                is_terminal_hydroxyl = True
                break
    if not is_terminal_hydroxyl:
        return False, "Hydroxyl group is not at terminal position"

    # Find the shortest path between carboxyl carbon and hydroxyl carbon
    path = Chem.rdmolops.GetShortestPath(mol, carboxyl_c_idx, hydroxyl_c_idx)
    if not path or len(path) == 0:
        return False, "No path between carboxyl and hydroxyl groups"

    # Check that the path is mostly carbons (allowing for oxygens)
    for idx in path:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() not in [6, 8]:  # Allow carbon and oxygen
            return False, "Path contains atoms other than carbon and oxygen"

    # Check that internal carbons are not branched
    for idx in path[1:-1]:  # Exclude terminal carbons
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 6:
            heavy_atom_neighbors = [a for a in atom.GetNeighbors() if a.GetAtomicNum() > 1]
            if len(heavy_atom_neighbors) > 2:
                return False, "Chain is branched between carboxyl and hydroxyl groups"

    return True, "Molecule is an omega-hydroxy fatty acid"

__metadata__ = {
    'chemical_class': {
        'name': 'omega-hydroxy fatty acid',
        'definition': 'Any member of the class of naturally-occurring straight-chain fatty acids n carbon atoms long with a carboxyl group at position 1 and a hydroxyl at position n (omega).',
    }
}