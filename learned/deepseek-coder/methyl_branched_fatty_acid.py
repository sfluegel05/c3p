"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
"""
Classifies: CHEBI:76971 methyl-branched fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl-branched fatty acid based on its SMILES string.
    A methyl-branched fatty acid is a branched-chain fatty acid containing methyl branches only.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl-branched fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-C(=O)O or -C(=O)[O-])
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0-,OX1H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for at least one methyl branch (C attached to exactly 3 other carbons)
    methyl_branch_pattern = Chem.MolFromSmarts("[CX4H3]")
    methyl_branch_matches = mol.GetSubstructMatches(methyl_branch_pattern)
    if len(methyl_branch_matches) == 0:
        return False, "No methyl branches found"

    # Check for long carbon chain (at least 4 carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 4:
        return False, "Carbon chain too short to be a fatty acid"

    # Ensure no non-methyl branches (e.g., ethyl, propyl, etc.)
    # We look for carbons with more than 2 connections that are not part of the main chain
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetDegree() > 2:
            # Check if the branch is longer than a single carbon (methyl group)
            neighbors = atom.GetNeighbors()
            for neighbor in neighbors:
                if neighbor.GetAtomicNum() == 6 and neighbor.GetDegree() > 1:
                    return False, "Non-methyl branches detected"

    return True, "Contains carboxylic acid group, methyl branches, and a long carbon chain"