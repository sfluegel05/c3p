"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid based on its SMILES string.
    A cyclic fatty acid is a fatty acid that contains at least one ring structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclic fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-C(=O)O or -COOH) or ester/amide
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    if not (mol.HasSubstructMatch(carboxylic_acid_pattern) or mol.HasSubstructMatch(ester_pattern)):
        return False, "No carboxylic acid or ester group found"

    # Check for the presence of rings
    ring_info = mol.GetRingInfo()
    if not ring_info.NumRings():
        return False, "No ring structure found"

    # Check for long carbon chain (at least 10 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
        return False, "Not enough carbons to be a fatty acid"

    # Check that the ring is part of the carbon chain
    # We can do this by ensuring that the ring atoms are connected to the carboxylic acid/ester group
    ring_atoms = set()
    for ring in ring_info.AtomRings():
        ring_atoms.update(ring)
    
    # Get atoms connected to the carboxylic acid/ester group
    acid_atoms = set()
    for match in mol.GetSubstructMatches(carboxylic_acid_pattern):
        acid_atoms.update(match)
    for match in mol.GetSubstructMatches(ester_pattern):
        acid_atoms.update(match)
    
    # Check if any ring atom is connected to the carboxylic acid/ester group
    if not ring_atoms.intersection(acid_atoms):
        return False, "Ring structure not connected to carboxylic acid/ester group"

    # If all conditions are met, it is a cyclic fatty acid
    return True, "Contains a carboxylic acid or ester group, a long carbon chain, and a ring structure connected to the acid/ester group"