"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
"""
Classifies: CHEBI:78231 N-acylphytosphingosine
"""
from rdkit import Chem

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine based on its SMILES string.
    An N-acylphytosphingosine is a ceramide with a phytosphingosine backbone (at least three hydroxyl groups)
    and a fatty acyl group (minimum 8 carbons) attached via an amide bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is N-acylphytosphingosine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Find amide groups [N connected to carbonyl]
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide group found"

    for amide_match in amide_matches:
        nitrogen_idx = amide_match[0]
        carbonyl_idx = amide_match[1]

        # Count acyl chain carbons (minimum 8)
        acyl_carbons = count_acyl_carbons(mol, carbonyl_idx, nitrogen_idx)
        if acyl_carbons < 8:
            continue

        # Count hydroxyls in sphingoid main chain (minimum 3)
        hydroxyls = get_hydroxyl_count_in_main_chain(mol, nitrogen_idx, carbonyl_idx)
        if hydroxyls >= 3:
            return True, f"Phytosphingosine backbone with {hydroxyls} hydroxyls and {acyl_carbons}-carbon acyl chain"

    return False, "Does not meet N-acylphytosphingosine criteria"

def count_acyl_carbons(mol, start_idx, exclude_idx):
    """Counts all carbons in the acyl chain connected to the carbonyl, excluding the nitrogen side."""
    visited = set()
    stack = [start_idx]
    carbon_count = 0

    while stack:
        atom_idx = stack.pop()
        if atom_idx in visited:
            continue
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() == 6:
            carbon_count += 1
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx == exclude_idx:
                continue
            if neighbor_idx not in visited:
                stack.append(neighbor_idx)
    return carbon_count

def get_hydroxyl_count_in_main_chain(mol, start_idx, exclude_idx):
    """Finds maximum hydroxyl count along the longest path in sphingoid base."""
    visited = set()
    max_hydroxyls = 0

    def dfs(current_idx, current_hydroxyls):
        nonlocal max_hydroxyls
        if current_idx in visited:
            return
        visited.add(current_idx)
        
        atom = mol.GetAtomWithIdx(current_idx)
        # Check for hydroxyl groups on this atom
        has_hydroxyl = False
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() >= 1:
                has_hydroxyl = True
                break
        
        new_hydroxyls = current_hydroxyls + (1 if has_hydroxyl else 0)
        if new_hydroxyls > max_hydroxyls:
            max_hydroxyls = new_hydroxyls
        
        # Explore all paths except back to exclude_idx
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx != exclude_idx and neighbor_idx not in visited:
                dfs(neighbor_idx, new_hydroxyls)
        
        visited.remove(current_idx)

    dfs(start_idx, 0)
    return max_hydroxyls