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

    # Check for presence of oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count == 0:
        return False, "No oxygen atoms found; not an oxygenated carotenoid"

    # Build a graph of conjugated atoms connected via conjugated bonds
    conj_atom_graph = {}  # atom_idx -> set of neighbor atom_idx
    for atom in mol.GetAtoms():
        if atom.GetIsConjugated():
            idx = atom.GetIdx()
            conj_atom_graph[idx] = set()

    # Build adjacency list for conjugated atoms
    for bond in mol.GetBonds():
        if bond.GetIsConjugated():
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetIsConjugated() and end_atom.GetIsConjugated():
                idx1 = begin_atom.GetIdx()
                idx2 = end_atom.GetIdx()
                conj_atom_graph[idx1].add(idx2)
                conj_atom_graph[idx2].add(idx1)

    # Find the largest conjugated system using connected components
    visited = set()
    max_conj_system_size = 0

    for atom_idx in conj_atom_graph:
        if atom_idx not in visited:
            # Perform DFS to find all connected conjugated atoms
            stack = [atom_idx]
            current_component_size = 0
            while stack:
                current_atom = stack.pop()
                if current_atom not in visited:
                    visited.add(current_atom)
                    current_component_size += 1
                    neighbors = conj_atom_graph[current_atom]
                    stack.extend(neighbors - visited)
            max_conj_system_size = max(max_conj_system_size, current_component_size)

    # Check if the largest conjugated system is long enough to be a carotenoid
    if max_conj_system_size < 18:
        return False, f"Insufficient conjugated system size ({max_conj_system_size}); not a carotenoid"

    # Check that the molecule has sufficient carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 35:
        return False, f"Too few carbon atoms ({c_count}); not a carotenoid"

    # Check for functional groups typical of xanthophylls
    has_oxygenated_functional_group = False
    # Hydroxyl group
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[OX2H]")):
        has_oxygenated_functional_group = True
    # Carbonyl group
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C=O")):
        has_oxygenated_functional_group = True
    # Epoxide ring
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[C]1-[O]-[C]1")):
        has_oxygenated_functional_group = True
    # Ether group
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C-O-C")):
        has_oxygenated_functional_group = True

    if not has_oxygenated_functional_group:
        return False, "No typical xanthophyll functional groups found"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.2f} Da); not a carotenoid"

    return True, "Molecule is an oxygenated carotenoid (xanthophyll)"

__metadata__ = {   
    'chemical_class': {   
        'id': 'CHEBI:27530',
        'name': 'xanthophyll',
        'definition': 'A subclass of carotenoids consisting of the oxygenated carotenes.',
    },
}