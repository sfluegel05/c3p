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
    if not mol:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (using a specific explicit pattern)
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)[O;H1]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"

    # Check for methyl branches
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            neighbor_carbons = [nbr for nbr in atom.GetNeighbors() if nbr.GetSymbol() == 'C']
            
            # Identify branches where carbon has more than two bonded carbon neighbors
            if len(neighbor_carbons) > 2:  
                # Ensure it's branching due to methyl groups
                methyl_neighbors = [nbr for nbr in neighbor_carbons if len([x for x in nbr.GetNeighbors() if x.GetSymbol() == 'C']) == 1]
                if len(methyl_neighbors) < (len(neighbor_carbons) - 2):
                    return False, "Non-methyl branching detected"

    return True, "Contains carboxylic acid group and only methyl branches"