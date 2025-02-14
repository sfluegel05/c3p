"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
"""
Classifies: trienoic fatty acid
Definition: Any polyunsaturated fatty acid that contains three double bonds.
"""
from rdkit import Chem

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    A trienoic fatty acid is a fatty acid with a terminal carboxylic acid group and at least three
    carbon-carbon double bonds in its longest unbranched aliphatic chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trienoic fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    carbox_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)

    if not carbox_matches:
        return False, "No terminal carboxylic acid group found"

    # Assume the fatty acid chain starts from the alpha carbon (adjacent to carboxylic group)
    carbox_carbons = [match[0] for match in carbox_matches]  # Carbonyl carbon atoms
    alpha_carbons = []
    for carbox_c in carbox_carbons:
        atom = mol.GetAtomWithIdx(carbox_c)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in carbox_carbons:
                alpha_carbons.append(neighbor.GetIdx())

    if not alpha_carbons:
        return False, "No alpha carbon found adjacent to carboxylic group"

    # Traverse the chain from the alpha carbon
    visited = set()
    chain_atoms = []
    num_double_bonds = 0

    def traverse_chain(atom_idx):
        nonlocal num_double_bonds
        stack = [(atom_idx, None)]  # (current_atom_idx, previous_atom_idx)
        while stack:
            current_idx, previous_idx = stack.pop()
            if current_idx in visited:
                return False  # Detected a cycle (branching)
            visited.add(current_idx)
            chain_atoms.append(current_idx)
            current_atom = mol.GetAtomWithIdx(current_idx)
            neighbors = [n for n in current_atom.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIdx() != previous_idx]
            # Check for branching
            if len(neighbors) > 1:
                return False  # Branching detected
            for bond in current_atom.GetBonds():
                if bond.GetOtherAtomIdx(current_idx) == previous_idx:
                    continue
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6:
                    num_double_bonds += 1
            if neighbors:
                neighbor_idx = neighbors[0].GetIdx()
                stack.append((neighbor_idx, current_idx))
        return True

    # Start traversal
    for alpha_c in alpha_carbons:
        visited.clear()
        chain_atoms.clear()
        num_double_bonds = 0
        if traverse_chain(alpha_c):
            break  # Successful traversal without branching
    else:
        return False, "Chain is branched, fatty acids are typically unbranched"

    # Count carbons in the chain
    c_count = len(chain_atoms) + 1  # Including the carboxylic carbon
    if c_count < 12:
        return False, f"Chain contains {c_count} carbon atoms, which is too short for a typical fatty acid"

    # Check the number of double bonds
    if num_double_bonds < 3:
        return False, f"Contains {num_double_bonds} carbon-carbon double bonds in the chain, needs at least 3"

    return True, "Molecule is a trienoic fatty acid with a terminal carboxylic acid group and at least three C=C double bonds in the unbranched aliphatic chain"