"""
Classifies: CHEBI:27325 xanthophyll
"""
"""
Classifies: xanthophyll
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    A xanthophyll is an oxygenated carotenoid, consisting of oxygenated carotenes.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a xanthophyll, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Exclude molecules containing heteroatoms other than oxygen
    allowed_atoms = {1, 6, 8}  # H, C, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, "Contains heteroatoms other than oxygen; not a xanthophyll"

    # Check for presence of oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count == 0:
        return False, "No oxygen atoms found; not an oxygenated carotenoid"

    # Check for carotenoid backbone (long conjugated polyene chain)
    # Create a graph of conjugated double bonds
    bonds = mol.GetBonds()
    conj_double_bonds = [bond for bond in bonds if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and bond.GetIsConjugated()]
    if not conj_double_bonds:
        return False, "No conjugated double bonds found; not a carotenoid"

    # Build adjacency list of conjugated double bonds
    bond_graph = {}
    for bond in conj_double_bonds:
        idx = bond.GetIdx()
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        connected_bonds = []
        for neighbor in begin_atom.GetBonds() + end_atom.GetBonds():
            if neighbor.GetIdx() == idx:
                continue
            if neighbor.GetBondType() == Chem.rdchem.BondType.DOUBLE and neighbor.GetIsConjugated():
                connected_bonds.append(neighbor.GetIdx())
        bond_graph[idx] = connected_bonds

    # Find the longest path in the bond graph
    def dfs(bond_idx, visited):
        visited.add(bond_idx)
        lengths = [1]
        for neighbor_idx in bond_graph[bond_idx]:
            if neighbor_idx not in visited:
                lengths.append(1 + dfs(neighbor_idx, visited.copy()))
        return max(lengths)

    max_conj_chain_length = 0
    for bond in conj_double_bonds:
        length = dfs(bond.GetIdx(), set())
        max_conj_chain_length = max(max_conj_chain_length, length)

    if max_conj_chain_length < 7:
        return False, f"Insufficient conjugated double bonds ({max_conj_chain_length}); not a carotenoid"

    # Check that the molecule has sufficient carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 25:
        return False, f"Too few carbon atoms ({c_count}); not a carotenoid"

    # Check for functional groups typical of xanthophylls
    has_hydroxyl = mol.HasSubstructMatch(Chem.MolFromSmarts("[OX2H]"))      # Hydroxyl group
    has_keto = mol.HasSubstructMatch(Chem.MolFromSmarts("C=O"))             # Carbonyl group
    has_epoxide = mol.HasSubstructMatch(Chem.MolFromSmarts("[C]1-[O]-[C]1"))  # Epoxide ring
    has_either = mol.HasSubstructMatch(Chem.MolFromSmarts("C-O-C"))         # Ether group

    if not (has_hydroxyl or has_keto or has_epoxide or has_either):
        return False, "No typical xanthophyll functional groups found"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecular weight too low ({mol_wt:.2f} Da); not a carotenoid"

    return True, "Molecule is an oxygenated carotenoid (xanthophyll)"

__metadata__ = {   
    'chemical_class': {   
        'id': 'CHEBI:27530',
        'name': 'xanthophyll',
        'definition': 'A subclass of carotenoids consisting of the oxygenated carotenes.',
    },
}