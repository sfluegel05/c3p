"""
Classifies: CHEBI:27325 xanthophyll
"""
"""
Classifies: Xanthophylls (oxygenated carotenes), a subclass of carotenoids.
Heuristic: Xanthophyll molecules are expected to have:
  - A long conjugated polyene chain (several consecutively conjugated double bonds),
  - A large carbon backbone (we require at least about 20 carbons),
  - At least one oxygenated functional group (e.g. hydroxyl, carbonyl, or epoxide).
In addition, we attempt to “measure” the conjugated system by building a graph of atoms
joined by non‐aromatic conjugated double bonds and computing its longest simple path.
Note: This is a heuristic; false positives/negatives may still occur.
"""

from rdkit import Chem

def longest_conjugated_chain(mol):
    """
    Returns the length (in number of atoms) of the longest chain composed of
    conjugated double bonds (non-aromatic). We build a graph in which an edge is added
    if a bond is a double bond, is conjugated, and is not aromatic.
    """
    # Build adjacency list only over bonds meeting the criteria.
    adj = {}
    for bond in mol.GetBonds():
        # Check that bond is double, conjugated, and not aromatic.
        if bond.GetBondTypeAsDouble() == 2.0 and bond.GetIsConjugated() and (not bond.GetIsAromatic()):
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            adj.setdefault(a1, set()).add(a2)
            adj.setdefault(a2, set()).add(a1)
            
    # Use DFS to compute longest simple path in this graph.
    longest = 0
    visited_global = set()
    
    def dfs(node, visited):
        nonlocal longest
        # path length in terms of atoms visited so far
        current_length = len(visited)
        if current_length > longest:
            longest = current_length
        # Explore neighbors not yet in visited.
        for nbr in adj.get(node, []):
            if nbr not in visited:
                dfs(nbr, visited | {nbr})
    
    for node in adj.keys():
        dfs(node, {node})
    return longest

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    Xanthophylls are oxygenated carotenoids; they typically contain a long conjugated
    polyene system, a large carbon backbone (even if not as large as C30, we require ~C20+),
    and one or more oxygenated functional groups (e.g. hydroxyl, carbonyl, or epoxide).

    Args:
        smiles (str): SMILES string for the molecule.

    Returns:
        bool: True if the molecule is classified as a xanthophyll, False otherwise.
        str: A reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon and oxygen atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if carbon_count < 20:
        return False, f"Too few carbons ({carbon_count}) for a typical xanthophyll structure"
    if oxygen_count < 1:
        return False, "No oxygen atoms present, so not a xanthophyll"
    
    # Check for oxygenated functional groups using SMARTS patterns
    hydroxyl = Chem.MolFromSmarts("[OX2H]")     # hydroxyl
    carbonyl = Chem.MolFromSmarts("[CX3]=O")      # carbonyl 
    epoxide  = Chem.MolFromSmarts("[#6]-1-[#8]-[#6]-1")  # epoxide ring
    if not (mol.HasSubstructMatch(hydroxyl) or mol.HasSubstructMatch(carbonyl) or mol.HasSubstructMatch(epoxide)):
        return False, "No oxygenated functional groups detected"
    
    # Check for a long conjugated polyene chain.
    # First, look for a specific pattern as a quick check.
    polyene_pattern = Chem.MolFromSmarts("C=C-C=C-C=C")  # Three alternating C=C bonds
    if mol.HasSubstructMatch(polyene_pattern):
        conjugation_ok = True
    else:
        # As a fallback, compute the longest chain of connected conjugated double bonds.
        chain_length = longest_conjugated_chain(mol)
        if chain_length < 8:
            return False, "No sufficiently long conjugated polyene chain detected (longest chain: {} atoms)".format(chain_length)
        conjugation_ok = True
    
    if conjugation_ok:
        return True, "Molecule has a long conjugated polyene chain, sufficient carbon backbone, and oxygen functionalities consistent with xanthophyll structure"
    else:
        return False, "Conjugated system criteria not met"

# Example usage:
# test_smiles = "OC1CC(C(=C(C1)C)/C=C/C(/C)=C/C=C/C(/C)=C/C=C/C=C(/C=C/C=C(/C=C/C1=C(C)C[C@@H](O)CC1(C)C)"  # simplified example
# result, reason = is_xanthophyll(test_smiles)
# print(result, reason)

# If the molecule fails any of these checks, we return False with the appropriate message.