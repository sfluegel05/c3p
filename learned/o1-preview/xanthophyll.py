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

    # Check for carotenoid backbone (conjugated polyene chain)
    # Approximate by counting longest conjugated chain of double bonds between carbons
    def get_longest_conj_chain(mol):
        longest = 0
        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and bond.GetIsConjugated():
                visited = set()
                queue = [(bond, 1)]
                while queue:
                    current_bond, length = queue.pop(0)
                    longest = max(longest, length)
                    for neighbor_bond in current_bond.GetBeginAtom().GetBonds() + current_bond.GetEndAtom().GetBonds():
                        if neighbor_bond.GetIdx() != current_bond.GetIdx() and neighbor_bond.GetIsConjugated() and neighbor_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            if neighbor_bond.GetIdx() not in visited:
                                visited.add(neighbor_bond.GetIdx())
                                queue.append((neighbor_bond, length + 1))
        return longest

    max_conj_double_bonds = get_longest_conj_chain(mol)
    if max_conj_double_bonds < 7:
        return False, f"Insufficient conjugated double bonds ({max_conj_double_bonds}); not a carotenoid"

    # Check that the molecule has sufficient carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 25:
        return False, f"Too few carbon atoms ({c_count}); not a carotenoid"

    # Check for functional groups typical of xanthophylls
    has_hydroxyl = mol.HasSubstructMatch(Chem.MolFromSmarts("[OX2H]"))      # Hydroxyl group
    has_keto = mol.HasSubstructMatch(Chem.MolFromSmarts("C=O"))             # Carbonyl group
    has_epoxide = mol.HasSubstructMatch(Chem.MolFromSmarts("[C]1-[O]-[C]1"))  # Epoxide ring

    if not (has_hydroxyl or has_keto or has_epoxide):
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