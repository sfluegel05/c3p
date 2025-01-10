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

    # Look for carboxylic acid group (use flexible pattern to account for terminal placement)
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)[O]")
    if mol.HasSubstructMatch(carboxylic_pattern) is False:
        return False, "No carboxylic acid group found"

    # Check for methyl-only branches
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':  # Check only carbon atoms
            non_hydrogen_neighbors = [neighbor for neighbor in atom.GetNeighbors() if neighbor.GetSymbol() != 'H']
            # Check if any carbon has 3 or fewer non-H bonds, potentially indicating a branch
            if len(non_hydrogen_neighbors) > 3:
                return False, "Non-methyl branching detected"
            # Specifically look for methyl group branches
            if len(non_hydrogen_neighbors) == 3:
                # Two out of 3 neighbors should be carbons with only one bond (CH3)
                methyl_count = sum(1 for neighbor in non_hydrogen_neighbors if neighbor.GetSymbol() == 'C' and neighbor.GetDegree() == 1)
                if methyl_count < 2:
                    return False, "Non-methyl branching detected"

    return True, "Contains carboxylic acid group and only methyl branches"