"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
"""
Classifies: branched-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid based on its SMILES string.
    A branched-chain fatty acid is any fatty acid in which the parent hydrocarbon chain has one or more alkyl substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1H]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Get carboxylic acid group matches
    matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    # Assume the first carboxyl group found is the main one
    carboxyl_c_idx = matches[0][0]  # Index of the carboxyl carbon

    # Find the carbon attached to the carboxyl carbon (alpha carbon)
    carboxyl_c_atom = mol.GetAtomWithIdx(carboxyl_c_idx)
    neighbors = carboxyl_c_atom.GetNeighbors()
    chain_carbon = None
    for neighbor in neighbors:
        if neighbor.GetAtomicNum() == 6:  # Carbon atom
            chain_carbon = neighbor
            break
    if chain_carbon is None:
        return False, "No carbon chain attached to carboxyl carbon"
    
    # Traverse the carbon chain to identify branching points
    visited_atoms = set()
    branching_points = set()
    
    def traverse_chain(atom, previous_atom_idx=None):
        visited_atoms.add(atom.GetIdx())
        # Count number of carbon neighbors excluding the previous atom
        carbon_neighbors = [n for n in atom.GetNeighbors() 
                            if n.GetAtomicNum() == 6 and n.GetIdx() != previous_atom_idx]
        if len(carbon_neighbors) > 1:
            branching_points.add(atom.GetIdx())
        for neighbor in carbon_neighbors:
            if neighbor.GetIdx() not in visited_atoms:
                traverse_chain(neighbor, atom.GetIdx())
                
    traverse_chain(chain_carbon, carboxyl_c_idx)
    
    # Count total carbons in the chain
    chain_length = len(visited_atoms)
    if chain_length < 4:
        return False, f"Chain too short for a fatty acid (length {chain_length})"
    
    if len(branching_points) == 0:
        return False, "No branching points found in the carbon chain"
    
    return True, "Contains a branched aliphatic chain attached to a carboxylic acid group"