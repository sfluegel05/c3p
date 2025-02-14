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
    from rdkit.Chem import rdMolOps

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
    if c_count < 30:
        return False, f"Too few carbons ({c_count}) to be a xanthophyll"

    # Check for presence of oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count == 0:
        return False, "No oxygen atoms found, not a xanthophyll"

    # Find the largest conjugated system
    def get_largest_conjugated_system_size(mol):
        emol = Chem.EditableMol(mol)
        bonds_to_remove = []
        for bond in mol.GetBonds():
            if not (bond.GetIsConjugated() and bond.GetBondType() != Chem.rdchem.BondType.AROMATIC):
                bonds_to_remove.append(bond.GetIdx())
        # Remove non-conjugated bonds
        for bond_idx in sorted(bonds_to_remove, reverse=True):
            bond = mol.GetBondWithIdx(bond_idx)
            emol.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        conj_mol = emol.GetMol()
        frags = Chem.GetMolFrags(conj_mol, asMols=True)
        if not frags:
            return 0
        # Get the size (number of atoms) of the largest conjugated system
        largest_frag_size = max(frag.GetNumAtoms() for frag in frags)
        return largest_frag_size

    largest_conj_system_size = get_largest_conjugated_system_size(mol)
    if largest_conj_system_size < 15:
        return False, f"Longest conjugated system is too small ({largest_conj_system_size} atoms) for a xanthophyll"

    # Check for typical oxygen-containing functional groups
    has_hydroxyl = mol.HasSubstructMatch(Chem.MolFromSmarts('[OX2H]'))
    has_ketone = mol.HasSubstructMatch(Chem.MolFromSmarts('C=O'))
    has_epoxide = mol.HasSubstructMatch(Chem.MolFromSmarts('C1OC1'))
    has_ether = mol.HasSubstructMatch(Chem.MolFromSmarts('COC'))
    has_carboxylic_acid = mol.HasSubstructMatch(Chem.MolFromSmarts('C(=O)[OH]'))
    has_ester = mol.HasSubstructMatch(Chem.MolFromSmarts('C(=O)O'))

    if not (has_hydroxyl or has_ketone or has_epoxide or has_ether or has_carboxylic_acid or has_ester):
        return False, "No typical oxygen-containing functional groups found in xanthophyll"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.2f} Da) for a typical xanthophyll"

    # Verify that the molecule has ring structures (many carotenoids have rings at the ends)
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 1:
        return False, f"Expected at least 1 ring in structure, found {num_rings}"

    # All checks passed
    return True, "Molecule is likely a xanthophyll (oxygenated carotenoid)"

# Note: The function checks for the following:
# - Valid SMILES string
# - Does not contain unusual elements (only C, H, O)
# - Sufficient number of carbons (minimum 30)
# - Presence of oxygen atoms
# - Size of the largest conjugated system (minimum 15 atoms)
# - Presence of typical functional groups in xanthophylls (hydroxyls, ketones, epoxides, ethers, carboxylic acids, esters)
# - Molecular weight consistent with xanthophylls
# - Presence of ring structures typical for carotenoids

# If all criteria are met, the molecule is classified as a xanthophyll.