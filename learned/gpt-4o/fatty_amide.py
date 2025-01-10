"""
Classifies: CHEBI:29348 fatty amide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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
    
    # Analyze carbon chain structure
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    longest_chain = 0
    
    for atom in carbon_atoms:
        visited = set()
        chain_length = _dfs_longest_chain(atom, visited)
        if chain_length > longest_chain:
            longest_chain = chain_length
    
    # Check for a minimum chain length of 12 carbon atoms
    if longest_chain < 12:
        return False, "Carbon chain too short to be a fatty acid-derived amide"
    
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
        if neighbor.GetAtomicNum() == 6:
            chain_length = max(chain_length, 1 + _dfs_longest_chain(neighbor, visited))
    
    visited.remove(atom.GetIdx())
    return chain_length