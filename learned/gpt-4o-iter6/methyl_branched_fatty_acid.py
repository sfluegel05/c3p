"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
from rdkit import Chem

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl-branched fatty acid based on its SMILES string.
    A methyl-branched fatty acid must have a terminal carboxylic acid group and only methyl branches.

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

    # Check for terminal carboxylic acid group
    # Here, we assume the carboxylic group is terminal (and might appear in several tautomeric forms)
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)OH")
    carboxylic_count = len(mol.GetSubstructMatch(carboxylic_pattern))
    if carboxylic_count != 1:
        return False, f"Expected 1 terminal carboxylic acid group, found {carboxylic_count}"
    
    # Traverse molecule to ensure methyl-only branches
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            neighbors = atom.GetNeighbors()
            # Check if current carbon is a branch point but not the terminal carboxylic group
            if atom.GetDegree() > 3:
                return False, "Non-methyl branches detected"
            # Ensure only CH3-type neighbors aside from main chain
            non_hydrogen_neighbors = sum(1 for neighbor in neighbors if neighbor.GetSymbol() != 'H')
            methyl_branches = sum(1 for neighbor in neighbors if neighbor.GetSymbol() == 'C' and neighbor.GetDegree() == 1)
            # There could be a situation where a branch isn't methyl or terminal pattern faults
            if non_hydrogen_neighbors > 1 and (non_hydrogen_neighbors - methyl_branches) > 1:
                return False, "Non-methyl branch detected"
                
    return True, "Contains only methyl branches with a terminal carboxylic acid group"