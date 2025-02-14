"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
"""
Classifies: CHEBI:24996 octadecadienoic acid
"""

from rdkit import Chem

def is_octadecadienoic_acid(smiles: str):
    """
    Determines if a molecule is an octadecadienoic acid based on its SMILES string.
    An octadecadienoic acid is any straight-chain C18 polyunsaturated fatty acid
    having two C=C double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octadecadienoic acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxylic_acid = Chem.MolFromSmarts("C(=O)[O;H1]")
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"

    # Ensure the molecule is acyclic (no rings)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings; not a straight-chain fatty acid"

    # Build adjacency list for carbon atoms
    adjacency = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            idx = atom.GetIdx()
            adjacency[idx] = []
    for bond in mol.GetBonds():
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
            adjacency[begin_atom.GetIdx()].append((end_atom.GetIdx(), bond))
            adjacency[end_atom.GetIdx()].append((begin_atom.GetIdx(), bond))

    # Function to perform depth-first search to find the longest carbon chain
    def dfs(node, visited, length, double_bonds):
        visited.add(node)
        max_length = length
        max_double_bonds = double_bonds
        for neighbor_idx, bond in adjacency.get(node, []):
            if neighbor_idx not in visited:
                bond_order = bond.GetBondType()
                bond_is_double = (bond_order == Chem.rdchem.BondType.DOUBLE)
                new_double_bonds = double_bonds + (1 if bond_is_double else 0)
                l, d = dfs(neighbor_idx, visited.copy(), length + 1, new_double_bonds)
                if l > max_length:
                    max_length = l
                    max_double_bonds = d
                elif l == max_length and d > max_double_bonds:
                    max_double_bonds = d
        return max_length, max_double_bonds

    max_chain_length = 0
    max_chain_double_bonds = 0
    # Iterate over all carbon atoms to find the longest chain
    for node in adjacency:
        length, double_bonds = dfs(node, set(), 1, 0)
        if length > max_chain_length:
            max_chain_length = length
            max_chain_double_bonds = double_bonds
        elif length == max_chain_length and double_bonds > max_chain_double_bonds:
            max_chain_double_bonds = double_bonds

    if max_chain_length < 18:
        return False, f"Longest carbon chain is {max_chain_length}, expected at least 18"
    if max_chain_double_bonds != 2:
        return False, f"Number of C=C double bonds along the chain is {max_chain_double_bonds}, expected 2"

    return True, "Molecule has a straight-chain C18 fatty acid backbone with two C=C double bonds"