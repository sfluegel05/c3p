"""
Classifies: CHEBI:27325 xanthophyll
"""
"""
Classifies: Xanthophylls (oxygenated carotenoids)
Heuristic:
  - Must have at least ~20 carbons.
  - Must contain at least one oxygen (and specifically at least one oxygenated functional group such as hydroxyl, methoxy, carbonyl, or epoxide).
  - Must have a long conjugated polyene system.
    We first quickly search for a pattern (three conjugated C=C bonds). If that fails, we compute the
    longest chain of connected (non‐aromatic, conjugated) double bonds using a DFS. We now require a minimum
    chain length of 12 atoms.
Note:
  Xanthophylls are oxygenated carotenoids. This is a heuristic that may not catch all edge cases.
"""

from rdkit import Chem

def longest_conjugated_chain(mol):
    """
    Returns the length (in number of atoms) of the longest chain composed exclusively of
    non‐aromatic conjugated double bonds. We build an undirected graph where an edge is added 
    if a bond is a double bond, conjugated, and not aromatic, and then use DFS to find the 
    longest simple path.
    """
    # Build adjacency list over bonds that are double, conjugated, and not aromatic.
    adj = {}
    for bond in mol.GetBonds():
        if (bond.GetBondTypeAsDouble() == 2.0 and bond.GetIsConjugated() 
              and (not bond.GetIsAromatic())):
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()
            adj.setdefault(i, set()).add(j)
            adj.setdefault(j, set()).add(i)
    
    longest = 0
    # DFS function which tracks visited nodes in the current path.
    def dfs(node, visited):
        nonlocal longest
        current_length = len(visited)
        if current_length > longest:
            longest = current_length
        for nbr in adj.get(node, []):
            if nbr not in visited:
                dfs(nbr, visited | {nbr})
                
    for node in adj.keys():
        dfs(node, {node})
    return longest

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    Xanthophylls are oxygenated carotenoids – they must have an extended conjugated
    polyene chain, a large carbon backbone, and one or more oxygenated functional groups.
    
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
    
    # Count carbons and oxygens.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if carbon_count < 20:
        return False, f"Too few carbons ({carbon_count}) for a typical xanthophyll structure"
    if oxygen_count < 1:
        return False, "No oxygen atoms present, so not a xanthophyll"
    
    # Define oxygenated functional groups via SMARTS:
    # hydroxyl, methoxy, carbonyl, and an epoxide pattern.
    hydroxyl = Chem.MolFromSmarts("[OX2H]")          # hydroxyl (–OH)
    methoxy  = Chem.MolFromSmarts("[OX2][CH3]")        # methoxy (–OCH3)
    carbonyl = Chem.MolFromSmarts("[CX3]=O")           # carbonyl (C=O)
    epoxide  = Chem.MolFromSmarts("[#6]-1-[#8]-[#6]-1") # three-membered epoxide ring
    
    if not (mol.HasSubstructMatch(hydroxyl) or 
            mol.HasSubstructMatch(methoxy) or 
            mol.HasSubstructMatch(carbonyl) or 
            mol.HasSubstructMatch(epoxide)):
        return False, "No oxygenated functional groups detected"
    
    # Check for the presence of a long conjugated polyene chain.
    # First, try a quick substructure search for 3 consecutively conjugated C=C bonds.
    polyene_pattern = Chem.MolFromSmarts("C=C-C=C-C=C")
    if mol.HasSubstructMatch(polyene_pattern):
        conjugation_ok = True
    else:
        chain_length = longest_conjugated_chain(mol)
        if chain_length < 12:
            return False, "No sufficiently long conjugated polyene chain detected (longest chain: {} atoms)".format(chain_length)
        conjugation_ok = True
    
    if conjugation_ok:
        return True, "Molecule has a long conjugated polyene chain, sufficient carbon backbone, " \
                     "and oxygenated functional groups consistent with xanthophyll structure"
    return False, "Criteria not met"
    
# Example usage:
# test_smiles = "OC1CC(C(=C(C1)C)/C=C/C(/C)=C/C=C/C(/C)=C/C=C/C=C(/C=C/C=C(/C=C/C1=C(C)C[C@@H](O)CC1(C)C"  # shortened/exemplary SMILES
# result, reason = is_xanthophyll(test_smiles)
# print(result, reason)