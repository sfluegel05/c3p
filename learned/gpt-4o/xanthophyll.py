"""
Classifies: CHEBI:27325 xanthophyll
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    A xanthophyll is an oxygenated carotenoid, characterized by a long chain
    of conjugated double bonds and one or more oxygen-containing functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a xanthophyll, False otherwise
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for an extensive conjugated polyene chain
    # Adjust pattern for detecting longer conjugated systems
    polyene_pattern = Chem.MolFromSmarts("([#6]=[#6])-c1ccccc1")
    if not mol.HasSubstructMatch(polyene_pattern):
        return False, "No extensive conjugated system characteristic of carotenoids found"

    # Check for common oxygen-containing groups in xanthophylls
    oxygen_patterns = [
        Chem.MolFromSmarts("[OX2H]"),  # Hydroxyl
        Chem.MolFromSmarts("[CX3]=[OX1]"),  # Carbonyl
        Chem.MolFromSmarts("O([CX4])[CX3]"),
        Chem.MolFromSmarts("[OX2]([#6])[#6]"),  # Ethers and more bonds
    ]
    if not any(mol.HasSubstructMatch(oxygen_pattern) for oxygen_pattern in oxygen_patterns):
        return False, "No matching oxygen-containing functional groups found"

    # Ensure sufficient number of double bonds
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE)
    atom_count = mol.GetNumAtoms()
    if double_bond_count < atom_count / 3:  # Proportional criterion for xanthophylls
        return False, f"Insufficient number of double bonds for carotenoid-like structure; found {double_bond_count}"

    return True, "Matches the structure of a xanthophyll with conjugated polyene and oxygen functionality"