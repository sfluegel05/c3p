"""
Classifies: CHEBI:134363 tertiary amine oxide
"""
"""
Classifies: tertiary amine N-oxide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_organic_group(mol, atom):
    """Helper function to check if an atom is part of an organic group"""
    # Consider atom organic if it's carbon or part of common organic substituents
    if atom.GetAtomicNum() == 6:  # Carbon
        return True
    if atom.GetAtomicNum() == 7:  # Nitrogen in organic group
        return True
    if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 0:  # Oxygen in ether/ester
        return True
    return False

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

    # Look for N+-O- pattern
    # Match any N+ connected to O- where N has 4 bonds total
    n_oxide_pattern = Chem.MolFromSmarts("[NX4+]-[O-]")
    if not mol.HasSubstructMatch(n_oxide_pattern):
        return False, "No N-oxide group found"
    
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
        neighbors = []
        for neighbor in n_atom.GetNeighbors():
            if neighbor.GetIdx() != o_idx:
                neighbors.append(neighbor)
        
        # Must have exactly 3 non-oxide neighbors
        if len(neighbors) != 3:
            continue
            
        # Check that all neighbors are part of organic groups
        organic_count = sum(1 for neighbor in neighbors if is_organic_group(mol, neighbor))
        if organic_count != 3:
            continue
            
        # Check that nitrogen has no hydrogens
        if n_atom.GetTotalNumHs() != 0:
            continue
            
        # Exclude N-oxides that are part of very small rings (likely not stable)
        ring_size_pattern = Chem.MolFromSmarts("[NX4+]1-[O-]-[*]1")
        if mol.HasSubstructMatch(ring_size_pattern):
            continue
            
        # All checks passed
        return True, "Found nitrogen with three organic groups and N-oxide"
    
    return False, "No valid tertiary amine N-oxide found"