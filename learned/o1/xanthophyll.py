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
    A xanthophyll is an oxygen-containing carotenoid; that is, a carotene (C40 polyene hydrocarbon) 
    with oxygen-containing functional groups such as hydroxyls, ketones, epoxides, etc.

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

    # Exclude molecules with unusual elements (only consider C, H, O)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (1, 6, 8):
            return False, f"Contains unusual element: {atom.GetSymbol()}"

    # Count number of carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, f"Too few carbons ({c_count}) to be a carotenoid"

    # Check for presence of oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count == 0:
        return False, "No oxygen atoms found, not a xanthophyll"

    # Find the longest chain of conjugated non-aromatic double bonds
    def get_longest_conjugated_chain(mol):
        longest_length = 0
        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and not bond.GetIsAromatic():
                length = dfs_conjugated_chain(bond, set())
                if length > longest_length:
                    longest_length = length
        return longest_length

    def dfs_conjugated_chain(bond, visited):
        bond_idx = bond.GetIdx()
        if bond_idx in visited or bond.GetIsAromatic():
            return 0
        visited.add(bond_idx)
        current_atom = bond.GetEndAtom()
        length = 1
        for nbr_bond in current_atom.GetBonds():
            nbr_bond_type = nbr_bond.GetBondType()
            if nbr_bond.GetIdx() != bond_idx and nbr_bond_type == Chem.rdchem.BondType.DOUBLE and not nbr_bond.GetIsAromatic():
                length += dfs_conjugated_chain(nbr_bond, visited)
                break  # Only consider linear conjugation
        return length

    longest_conj_chain = get_longest_conjugated_chain(mol)
    if longest_conj_chain < 7:
        return False, f"Longest conjugated chain is too short ({longest_conj_chain}) for a carotenoid"

    # Check for typical oxygen-containing functional groups
    has_hydroxyl = mol.HasSubstructMatch(Chem.MolFromSmarts('[OX2H]'))
    has_ketone = mol.HasSubstructMatch(Chem.MolFromSmarts('C=O'))
    has_epoxide = mol.HasSubstructMatch(Chem.MolFromSmarts('C1OC1'))
    has_ether = mol.HasSubstructMatch(Chem.MolFromSmarts('COC'))

    if not (has_hydroxyl or has_ketone or has_epoxide or has_ether):
        return False, "No typical oxygen-containing functional groups found in xanthophyll"

    # Check for carotenoid-like structure: long linear or semi-linear molecule
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.2f} Da) for a typical xanthophyll"

    # Verify that the molecule is not cyclic (many carotenoids have ring structures at the ends)
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 2:
        return False, f"Expected at least 2 rings in structure, found {num_rings}"

    # All checks passed
    return True, "Molecule is likely a xanthophyll (oxygenated carotenoid)"

# Note: The function checks for the following:
# - Valid SMILES string
# - Does not contain unusual elements (only C, H, O)
# - Sufficient number of carbons (minimum 20)
# - Presence of oxygen atoms
# - Longest chain of conjugated non-aromatic double bonds (minimum length 7)
# - Presence of typical functional groups in xanthophylls (hydroxyls, ketones, epoxides, ethers)
# - Molecular weight consistent with xanthophylls
# - Presence of ring structures typical for carotenoids

# If all criteria are met, the molecule is classified as a xanthophyll.