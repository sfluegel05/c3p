"""
Classifies: CHEBI:29348 fatty amide
"""
from rdkit import Chem

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide based on its SMILES string.
    A fatty amide is a monocarboxylic acid amide derived from a fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty amide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the amide group pattern
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide group found"

    # Identify all carbon chains - flexible to account for both straight and branched
    carbon_chains = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            chain_length = _dfs_longest_chain(atom, set())
            if chain_length > 0:
                carbon_chains.append(chain_length)

    # Check for at least one long carbon chain typical of fatty acids
    if not any(chain >= 12 for chain in carbon_chains):
        return False, "No sufficiently long carbon chain detected"
    
    return True, "Contains an amide group and a long carbon chain characteristic of fatty amides"

def _dfs_longest_chain(atom, visited):
    """
    Depth-First Search to find the longest chain of connected carbon atoms.

    Args:
        atom (rdkit.Chem.rdchem.Atom): The starting atom for the DFS
        visited (set): The set of visited atoms to avoid cycles

    Returns:
        int: Length of the longest chain found starting from the given atom
    """
    if atom.GetIdx() in visited:
        return 0
    visited.add(atom.GetIdx())
    
    chain_length = 1  # Count the current atom
    for neighbor in atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:  # Check only carbon neighbors for the chain
            chain_length = max(chain_length, 1 + _dfs_longest_chain(neighbor, visited))
    
    visited.remove(atom.GetIdx())
    return chain_length