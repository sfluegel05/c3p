"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
"""
Classifies: omega-hydroxy fatty acid
"""

from rdkit import Chem

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

    # Identify carboxyl group [C(=O)[O-,OH]]
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1H0-,OX2H1]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No carboxyl group found"

    # Identify hydroxyl group [O-H or O-]
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H,OX1-]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if not hydroxyl_matches:
        return False, "No hydroxyl group found"

    # Loop over all carboxyl and hydroxyl groups to find a matching pair
    for carb_match in carboxyl_matches:
        carboxyl_c_idx = carb_match[0]  # Carbon of carboxyl group
        carboxyl_c_atom = mol.GetAtomWithIdx(carboxyl_c_idx)

        # Check if carboxyl carbon is connected to a carbon chain
        chain_start_atoms = [nbr for nbr in carboxyl_c_atom.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != carboxyl_c_idx]
        if not chain_start_atoms:
            continue  # Carboxyl group is not connected to a carbon chain
        chain_start_idx = chain_start_atoms[0].GetIdx()

        for hydro_match in hydroxyl_matches:
            hydroxyl_o_idx = hydro_match[0]
            hydroxyl_o_atom = mol.GetAtomWithIdx(hydroxyl_o_idx)

            # Identify the carbon attached to the hydroxyl group
            hydroxyl_c_atoms = [nbr for nbr in hydroxyl_o_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
            if not hydroxyl_c_atoms:
                continue  # Hydroxyl group is not attached to a carbon
            hydroxyl_c_atom = hydroxyl_c_atoms[0]
            hydroxyl_c_idx = hydroxyl_c_atom.GetIdx()

            # Ensure the hydroxyl carbon is terminal (connected to only one heavy atom besides the oxygen)
            hydroxyl_c_neighbors = [nbr for nbr in hydroxyl_c_atom.GetNeighbors() if nbr.GetAtomicNum() > 1 and nbr.GetIdx() != hydroxyl_o_idx]
            if len(hydroxyl_c_neighbors) > 1:
                continue  # Hydroxyl carbon is not terminal

            # Find the path between the carboxyl chain start and the hydroxyl carbon
            path = Chem.rdmolops.GetShortestPath(mol, chain_start_idx, hydroxyl_c_idx)
            if not path:
                continue  # No path found

            # Check that the path is linear and unbranched
            branched = False
            for idx in path:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() != 6:
                    branched = True  # Non-carbon atom in chain
                    break
                num_heavy_neighbors = sum(1 for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1)
                if (idx == chain_start_idx or idx == hydroxyl_c_idx) and num_heavy_neighbors > 2:
                    branched = True  # Terminal carbons should not be branched
                    break
                elif idx != chain_start_idx and idx != hydroxyl_c_idx and num_heavy_neighbors > 2:
                    branched = True  # Internal carbons should have max 2 heavy neighbors
                    break
            if branched:
                continue  # Not a linear chain

            # Check for rings in the path
            if mol.GetRingInfo().IsAtomInRingOfSize(idx, 3) or mol.GetRingInfo().IsAtomInRingOfSize(idx, 4):
                continue  # Ring structures are not allowed

            # All checks passed; molecule is an omega-hydroxy fatty acid
            return True, "Molecule is an omega-hydroxy fatty acid"

    return False, "Molecule does not match criteria for an omega-hydroxy fatty acid"

__metadata__ = {
    'chemical_class': {
        'name': 'omega-hydroxy fatty acid',
        'definition': 'Any member of the class of naturally-occurring straight-chain fatty acids n carbon atoms long with a carboxyl group at position 1 and a hydroxyl at position n (omega).',
    }
}