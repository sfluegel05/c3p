"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    An aliphatic nitrile is a nitrile derived from an aliphatic compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic nitrile, False otherwise.
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for nitrile group (C#N)
    nitrile_pattern = Chem.MolFromSmarts("C#N")
    if not mol.HasSubstructMatch(nitrile_pattern):
        return False, "No nitrile group found"

    # Check for the presence of metal atoms
    for atom in mol.GetAtoms():
      if atom.GetAtomicNum() not in [1, 6, 7, 8, 9, 15, 16, 17, 35, 53]:
        return False, "Contains a non-organic atom"

    # Check for aromatic rings
    num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if num_aromatic_rings > 0:
        return False, "Aromatic ring detected, not an aliphatic nitrile."

    return True, "Contains a nitrile group, no aromatic rings, and only allowed elements."