"""
Classifies: CHEBI:27325 xanthophyll
"""
"""
Classifies: Xanthophylls (oxygenated carotenoids)
Heuristic:
  - Must have a large carbon backbone (at least 30 carbons) – typical carotenoids have ~40.
  - Must have at least one oxygenated functional group.
      We first check common groups (hydroxyl, methoxy, carbonyl, epoxide) and if those fail,
      we check if any oxygen is attached to an sp² (conjugated) carbon.
  - Must contain a long, conjugated polyene chain.
    We first look for a quick pattern of three consecutive C=C bonds; if that fails,
    we compute the longest chain of connected conjugated double bonds (non‐aromatic) using DFS
    and require at least 12 atoms in the chain.
Note:
  This is an heuristic to capture oxygenated carotenoids (xanthophylls) and may have false positives/negatives.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_conjugated_chain(mol):
    """
    Returns the length (number of atoms) of the longest chain composed exclusively of
    non‐aromatic, conjugated double bonds. We construct an undirected graph where an edge is added 
    if a bond is a double bond, conjugated, and not aromatic, and then use DFS to find the longest simple path.
    """
    # Build adjacency list over bonds that are double, conjugated, and not aromatic.
    adj = {}
    for bond in mol.GetBonds():
        # Check if the bond is a double bond, conjugated, and not aromatic.
        if bond.GetBondTypeAsDouble() == 2.0 and bond.GetIsConjugated() and not bond.GetIsAromatic():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()
            adj.setdefault(i, set()).add(j)
            adj.setdefault(j, set()).add(i)
    
    longest = 0
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
    Xanthophylls are oxygenated carotenoids. They must:
      - Have a large number of carbons (≥30).
      - Contain one or more oxygenated functional groups. We check for common groups
        (hydroxyl, methoxy, carbonyl, epoxide) and if none are found, we check if any oxygen
        is attached to an sp²-hybridized (conjugated) carbon.
      - Possess a long conjugated polyene chain (quick match for three conjugated C=C bonds or a 
        computed chain of at least 12 atoms).
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if the molecule is classified as a xanthophyll.
      str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count the number of carbon atoms.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 30:
        return False, f"Too few carbons ({carbon_count}); expected at least 30 for a carotenoid skeleton"
    
    # Check if there is at least one oxygen atom.
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 1:
        return False, "No oxygen atoms present, so not a xanthophyll"
    
    # Check for common oxygenated functional groups using SMARTS.
    hydroxyl = Chem.MolFromSmarts("[OX2H]")          # hydroxyl (–OH)
    methoxy  = Chem.MolFromSmarts("[OX2][CH3]")        # methoxy (–OCH3)
    carbonyl = Chem.MolFromSmarts("[CX3]=O")           # carbonyl (C=O)
    epoxide  = Chem.MolFromSmarts("[#6]-1-[#8]-[#6]-1") # epoxide (three-membered ring)
    
    has_oxygen_group = (mol.HasSubstructMatch(hydroxyl) or 
                        mol.HasSubstructMatch(methoxy) or 
                        mol.HasSubstructMatch(carbonyl) or 
                        mol.HasSubstructMatch(epoxide))
    
    # If no common group was found, try a looser check:
    if not has_oxygen_group:
        # Check if any oxygen is attached directly to an sp2 (conjugated) carbon.
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 8:  # oxygen atom
                for nbr in atom.GetNeighbors():
                    if nbr.GetAtomicNum() == 6 and nbr.GetHybridization() in [Chem.rdchem.HybridizationType.SP2, Chem.rdchem.HybridizationType.SP]:
                        # Instead of calling GetIsConjugated on the atom, check the connecting bond.
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                        if bond and bond.GetIsConjugated():
                            has_oxygen_group = True
                            break
                if has_oxygen_group:
                    break
    if not has_oxygen_group:
        return False, "No oxygenated functional groups detected"
    
    # Check for a long conjugated polyene chain.
    # First, try a quick SMARTS search for three consecutive conjugated C=C bonds.
    polyene_pattern = Chem.MolFromSmarts("C=C-C=C-C=C")
    if mol.HasSubstructMatch(polyene_pattern):
        conjugation_ok = True
    else:
        chain_length = longest_conjugated_chain(mol)
        if chain_length < 12:
            return False, f"No sufficiently long conjugated polyene chain detected (longest chain: {chain_length} atoms)"
        conjugation_ok = True

    if conjugation_ok:
        return True, "Molecule has a long conjugated polyene chain, a large carbon skeleton (≥30 carbons), and oxygenated groups consistent with a xanthophyll structure"
    return False, "Criteria not met"

# Example usage:
# test_smiles = "CC(/C=C/C=C(C)/C=C/C1=C(C)CCCC1(C)C)=C\\C=C\\C=C(C)\\C=C\\C=C(/C)C1OC2(C)CCCC(C)(C)C2=C1"
# result, reason = is_xanthophyll(test_smiles)
# print(result, reason)