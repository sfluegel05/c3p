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

    # Check for carboxylic acid group using a specific SMARTS pattern
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)[O;H1]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"

    # Check for presence of methyl branches
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            # Neighboring carbons
            neighbor_carbons = [nbr for nbr in atom.GetNeighbors() if nbr.GetSymbol() == 'C']
            
            # Look for carbon atoms that should act as branch points
            branch_neighbors = [nbr for nbr in neighbor_carbons if len(nbr.GetNeighbors()) > 1]

            # More than two neighbors at a single carbon indicates a potential branch
            if len(neighbor_carbons) > 2:
                # Ensure that all branch neighbors are methyl groups
                for branch_neighbor in branch_neighbors:
                    # Methyl group should have 1 carbon and 3 hydrogens
                    if sum(1 for x in branch_neighbor.GetNeighbors() if x.GetSymbol() in ['C', 'H']) != 4:
                        return False, "Non-methyl branching detected"

    return True, "Contains carboxylic acid group and only methyl branches"