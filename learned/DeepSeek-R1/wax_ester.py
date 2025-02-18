"""
Classifies: CHEBI:10036 wax ester
"""
"""
Classifies: CHEBI:100157 wax ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a wax ester based on its SMILES string.
    A wax ester is an ester of a fatty acid and a fatty alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # Must have exactly one ester group
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"
    
    # Get the oxygen and carbonyl carbon from the ester group
    ester_atoms = ester_matches[0]
    oxygen_idx = ester_atoms[0]
    carbonyl_carbon_idx = ester_atoms[1]
    
    # Get neighboring atoms to determine chains
    oxygen = mol.GetAtomWithIdx(oxygen_idx)
    carbonyl_carbon = mol.GetAtomWithIdx(carbonyl_carbon_idx)
    
    # Alcohol chain: connected to oxygen (excluding the ester)
    alcohol_neighbors = [n for n in oxygen.GetNeighbors() if n.GetIdx() != carbonyl_carbon_idx]
    if not alcohol_neighbors:
        return False, "No alcohol chain attached to ester oxygen"
    alcohol_start = alcohol_neighbors[0]
    
    # Acid chain: connected to carbonyl carbon (excluding the oxygen)
    acid_neighbors = [n for n in carbonyl_carbon.GetNeighbors() if n.GetIdx() != oxygen_idx]
    if not acid_neighbors:
        return False, "No acid chain attached to ester carbonyl"
    acid_start = acid_neighbors[0]
    
    # Helper function to count carbons in a chain starting from atom, excluding ester atoms
    def count_chain_carbons(atom, visited=None):
        if visited is None:
            visited = set()
        if atom.GetIdx() in visited:
            return 0
        if atom.GetAtomicNum() != 6:  # Only consider carbons
            return 0
        visited.add(atom.GetIdx())
        count = 1  # Count this carbon
        for neighbor in atom.GetNeighbors():
            # Avoid revisiting ester atoms (oxygen and carbonyl carbon)
            if neighbor.GetIdx() in {oxygen_idx, carbonyl_carbon_idx}:
                continue
            count += count_chain_carbons(neighbor, visited)
        return count
    
    # Count carbons in alcohol and acid chains
    alcohol_carbons = count_chain_carbons(alcohol_start)
    acid_carbons = count_chain_carbons(acid_start)
    
    # Check chain lengths (minimum 8 carbons each for typical wax esters)
    if alcohol_carbons < 8 or acid_carbons < 8:
        return False, f"Chains too short (alcohol: {alcohol_carbons}, acid: {acid_carbons})"
    
    return True, "Contains a single ester group linking two long hydrocarbon chains"