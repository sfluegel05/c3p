"""
Classifies: CHEBI:26208 polyunsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_polyunsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a polyunsaturated fatty acid based on its SMILES string.
    A polyunsaturated fatty acid is a carboxylic acid with a carbon chain containing more than one
    aliphatic (non-aromatic) double bond. The chain may have minor branching or functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyunsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify terminal carboxylic acid group (COOH)
    carboxylic_acid = Chem.MolFromSmarts('C(=O)[O;H1]')
    carb_matches = mol.GetSubstructMatches(carboxylic_acid)
    if not carb_matches:
        return False, "No carboxylic acid group found"

    # Assume the first match is the carboxylic acid group
    carb_atom_idx = carb_matches[0][0]  # Index of the carbonyl carbon

    # Check if the carboxyl carbon is connected to a carbon chain (aliphatic chain)
    carboxyl_carbon = mol.GetAtomWithIdx(carb_atom_idx)
    connected_carbons = [neighbor for neighbor in carboxyl_carbon.GetNeighbors() if neighbor.GetAtomicNum() == 6]
    if not connected_carbons:
        return False, "Carboxyl carbon is not connected to carbon chain"

    # Use BFS to traverse the entire carbon chain connected to the carboxyl carbon
    visited = set()
    queue = []
    double_bond_count = 0
    carbon_chain_length = 0
    aromatic_atoms = set()
    
    # Initialize BFS with carbons connected to carboxyl carbon
    for neighbor in connected_carbons:
        if neighbor.GetAtomicNum() == 6:
            queue.append((neighbor, carboxyl_carbon.GetIdx()))
            visited.add(neighbor.GetIdx())

    while queue:
        current_atom, previous_atom_idx = queue.pop(0)
        carbon_chain_length += 1  # Count the carbons in the chain

        # Check if the atom is in an aromatic ring
        if current_atom.GetIsAromatic():
            aromatic_atoms.add(current_atom.GetIdx())

        # Iterate over neighbors
        for bond in current_atom.GetBonds():
            neighbor = bond.GetOtherAtom(current_atom)
            neighbor_idx = neighbor.GetIdx()

            # Ignore hydrogen atoms
            if neighbor.GetAtomicNum() != 6:
                continue

            # Avoid going back to the previous atom
            if neighbor_idx == previous_atom_idx:
                continue

            # If we have already visited this atom, skip it
            if neighbor_idx in visited:
                continue

            # Add neighbor to queue
            queue.append((neighbor, current_atom.GetIdx()))
            visited.add(neighbor_idx)

            # Check if the bond is a double bond and non-aromatic
            if bond.GetBondType() == rdchem.BondType.DOUBLE and not bond.GetIsAromatic():
                double_bond_count += 1

    # Exclude molecules with aromatic rings
    if aromatic_atoms:
        return False, "Molecule contains aromatic rings"

    # Check if there are more than one double bond
    if double_bond_count > 1:
        return True, f"Contains {double_bond_count} aliphatic double bonds in carbon chain"
    else:
        return False, f"Contains {double_bond_count} aliphatic double bond(s); requires more than one"