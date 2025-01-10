"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl-branched fatty acid based on its SMILES string.
    A methyl-branched fatty acid must have a carboxylic acid group and only methyl branches.

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

    # Look for the carboxylic acid group (C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Look for methyl branches (CH3-)
    # A methyl group is a carbon with three hydrogens
    methyl_group_pattern = Chem.MolFromSmarts("[CH3]")
    methyl_matches = mol.GetSubstructMatches(methyl_group_pattern)
    if not methyl_matches:
        return False, "No methyl branches detected"

    # Check if all branches from the main chain are methyl
    for atom in mol.GetAtoms():
        # If atom is carbon and has branches
        if atom.GetSymbol() == 'C' and atom.GetDegree() > 2:
            num_methyl_branch = 0
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C' and neighbor.GetDegree() == 1:
                    num_methyl_branch += 1
            # If any non-methyl branches detected, return False
            if (atom.GetDegree() - num_methyl_branch) > 2:
                return False, "Non-methyl branches detected"

    return True, "Contains methyl branches only with a carboxylic acid group"