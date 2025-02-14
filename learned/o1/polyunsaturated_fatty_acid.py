"""
Classifies: CHEBI:26208 polyunsaturated fatty acid
"""
from rdkit import Chem

def is_polyunsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a polyunsaturated fatty acid based on its SMILES string.
    A polyunsaturated fatty acid is a long-chain aliphatic carboxylic acid containing more than one double bond.

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
    carboxylic_acid = Chem.MolFromSmarts("C(=O)[O;H1]")
    carb_matches = mol.GetSubstructMatches(carboxylic_acid)
    if not carb_matches:
        return False, "No terminal carboxylic acid group found"

    # Assume the first match is the terminal carboxylic acid
    carb_atom_idx = carb_matches[0][0]  # Index of the carbonyl carbon

    # Use breadth-first search (BFS) to traverse the aliphatic chain
    from collections import deque

    visited = set()
    queue = deque()
    chain_atoms = set()

    # Start from the carbonyl carbon and explore connected carbons
    for neighbor in mol.GetAtomWithIdx(carb_atom_idx).GetNeighbors():
        if neighbor.GetAtomicNum() == 6 and not neighbor.IsInRing():
            queue.append(neighbor.GetIdx())
            visited.add(neighbor.GetIdx())

    while queue:
        atom_idx = queue.popleft()
        atom = mol.GetAtomWithIdx(atom_idx)
        chain_atoms.add(atom_idx)
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            # Continue if neighbor is carbon, not in a ring, and not visited
            if neighbor.GetAtomicNum() == 6 and not neighbor.IsInRing() and n_idx not in visited:
                queue.append(n_idx)
                visited.add(n_idx)

    # Check if the chain length is sufficient (e.g., at least 8 carbons)
    if len(chain_atoms) < 8:
        return False, f"Aliphatic chain too short ({len(chain_atoms)} carbons)"

    # Count carbon-carbon double bonds in the aliphatic chain
    double_bond_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            # Check if both atoms are carbons in the chain
            if (begin_idx in chain_atoms and end_idx in chain_atoms and
                mol.GetAtomWithIdx(begin_idx).GetAtomicNum() == 6 and
                mol.GetAtomWithIdx(end_idx).GetAtomicNum() == 6):
                # Exclude aromatic and ring bonds
                if not bond.GetIsAromatic() and not bond.IsInRing():
                    double_bond_count += 1

    if double_bond_count > 1:
        return True, f"Contains {double_bond_count} carbon-carbon double bonds in the aliphatic chain"
    else:
        return False, f"Contains {double_bond_count} double bond(s); requires more than one"