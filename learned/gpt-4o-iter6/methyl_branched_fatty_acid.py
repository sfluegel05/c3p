"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
from rdkit import Chem

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

    # Check for the carboxylic acid group (must be terminal, C(=O)O must not have further connections)
    terminal_carboxylic_pattern = Chem.MolFromSmarts("C(=O)[O;H1]")
    if not mol.HasSubstructMatch(terminal_carboxylic_pattern):
        return False, "No terminal carboxylic acid group found"
    
    # Look for non-methyl branching
    for atom in mol.GetAtoms():
        # If a carbon atom has more than one non-hydrogen neighbor, it could be a branching point
        if atom.GetSymbol() == 'C' and atom.GetDegree() > 3:
            # Count how many of the neighbors are methyl-type carbons (CH3)
            methyl_branches = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetSymbol() == 'C' and neighbor.GetDegree() == 1)
            # If there is a branch that is not a methyl group, return False
            if atom.GetDegree() - methyl_branches > 2:
                return False, "Non-methyl branches detected"

    return True, "Contains only methyl branches with a terminal carboxylic acid group"