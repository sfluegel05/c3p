"""
Classifies: CHEBI:134363 tertiary amine oxide
"""
"""
Classifies: tertiary amine N-oxide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_organic_group(mol, atom, visited=None):
    """
    Helper function to check if an atom is part of a valid organic group
    Uses recursive checking to validate the whole connected group
    """
    if visited is None:
        visited = set()
    
    # Prevent infinite recursion
    if atom.GetIdx() in visited:
        return True
    visited.add(atom.GetIdx())
    
    # Basic organic atoms
    if atom.GetAtomicNum() not in [6, 7, 8, 16]:  # C, N, O, S
        return False
        
    # Check neighbors recursively
    for neighbor in atom.GetNeighbors():
        if neighbor.GetIdx() not in visited:
            if not is_organic_group(mol, neighbor, visited):
                return False
                
    return True

def is_tertiary_amine_oxide(smiles: str):
    """
    Determines if a molecule is a tertiary amine N-oxide based on its SMILES string.
    A tertiary amine N-oxide has a nitrogen atom bonded to three organic groups and 
    an oxygen atom with a formal charge separation (N+-O-).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine N-oxide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for N+-O- pattern where N is not in a ring with O
    # and has exactly 4 bonds (3 to C + 1 to O)
    n_oxide_pattern = Chem.MolFromSmarts("[NX4+;!R]-[O-]")
    if not mol.HasSubstructMatch(n_oxide_pattern):
        return False, "No acyclic N-oxide group found"
    
    # Get matches for N-oxide groups
    n_oxide_matches = mol.GetSubstructMatches(n_oxide_pattern)
    
    for match in n_oxide_matches:
        n_idx = match[0]  # Index of nitrogen atom
        o_idx = match[1]  # Index of oxygen atom
        n_atom = mol.GetAtomWithIdx(n_idx)
        o_atom = mol.GetAtomWithIdx(o_idx)
        
        # Verify charges
        if n_atom.GetFormalCharge() != 1 or o_atom.GetFormalCharge() != -1:
            continue
            
        # Count non-oxide neighbors of nitrogen
        organic_neighbors = []
        for neighbor in n_atom.GetNeighbors():
            if neighbor.GetIdx() != o_idx:
                # Check if neighbor is start of a valid organic group
                if is_organic_group(mol, neighbor):
                    organic_neighbors.append(neighbor)
        
        # Must have exactly 3 organic neighbors
        if len(organic_neighbors) != 3:
            continue
            
        # Check that nitrogen has no hydrogens
        if n_atom.GetTotalNumHs() != 0:
            continue
            
        # Check that nitrogen is not part of a ring with the oxide oxygen
        ring_info = mol.GetRingInfo()
        if ring_info.AreAtomsBonded(n_idx, o_idx):
            continue
            
        # Check that the N-oxide is not part of a more complex functional group
        # like N-oxide heterocycles or N-oxide derivatives
        complex_pattern = Chem.MolFromSmarts("[NX4+](-[O-])-[!C]")
        if mol.HasSubstructMatch(complex_pattern):
            continue
            
        # All checks passed
        return True, "Found nitrogen with three organic groups and N-oxide"
    
    return False, "No valid tertiary amine N-oxide found"