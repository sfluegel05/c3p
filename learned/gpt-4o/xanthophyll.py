"""
Classifies: CHEBI:27325 xanthophyll
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    A xanthophyll is an oxygenated carotenoid, usually with a long chain of conjugated double bonds
    and one or more oxygen-containing functional groups.

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

    # Check for a conjugated polyene chain typical of carotenoids
    polyene_pattern = Chem.MolFromSmarts("[#6]=[#6]([#6])-[#6]=[#6]([#6])-[#6]=[#6]([#6])")
    if not mol.HasSubstructMatch(polyene_pattern):
        return False, "No conjugated polyene chain found"

    # Check for oxygen-containing functional groups, such as hydroxyl, carbonyl or epoxide
    oxy_patterns = [
        Chem.MolFromSmarts("[OX2H]"),  # Hydroxyl group
        Chem.MolFromSmarts("[CX3]=[OX1]"),  # Carbonyl group
        Chem.MolFromSmarts("[OX2][CX3]1[#6][#6][O][#6]1"),  # Epoxide group
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in oxy_patterns):
        return False, "No oxygen-containing functional groups found"

    # Optional check: Count the total number of double bonds to ensure a carotenoid-like structure
    num_double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE)
    if num_double_bonds < 4:  # Arbitrary threshold, many xanthophylls have more than 10
        return False, "Insufficient number of double bonds for carotenoid structure"

    return True, "Contains a conjugated polyene chain and oxygen-containing functional groups, indicative of xanthophyll"