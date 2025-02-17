"""
Classifies: CHEBI:67194 cannabinoid
"""
"""
Classifies: Cannabinoid chemical entities

Definition:
  "A diverse group of pharmacologically active secondary metabolites characteristic
   to the Cannabis plant as well as produced naturally in the body. Cannabinoids
   contain oxygen as part of a heterocyclic ring or in the form of various functional groups."

Heuristic rules in this version:
  1. Preliminary validity checks: molecule must be valid, contain oxygen, have at least 15 carbons
     and a molecular weight above ~200 Da.
  2. For phytocannabinoids: require a resorcinol-like motif (two SMARTS variants) with at least one alkyl side chain
     of moderate length (between 5 and 10 carbons) attached to an aromatic carbon.
  3. For synthetic cannabinoids: check for an indole ring using an improved SMARTS pattern.
  4. For endocannabinoids: if the molecule has a polar group (amide, ester or acid)
     and the longest overall carbon chain is at least 16 carbons then classify as cannabinoid.
  5. Fallback: if there is an oxygenated heterocycle (ring size ≥5 with at least one oxygen)
     and a long aliphatic chain (≥16 carbons), then flag as potentially cannabinoid.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_carbon_chain(mol):
    """
    Returns the length (number of atoms) of the longest continuous chain of carbon atoms.
    We build a connectivity graph for carbon only and then use DFS.
    """
    carbon_idxs = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    graph = {idx: [] for idx in carbon_idxs}
    for bond in mol.GetBonds():
        a = bond.GetBeginAtomIdx()
        b = bond.GetEndAtomIdx()
        if a in graph and b in graph:
            graph[a].append(b)
            graph[b].append(a)
            
    longest = 0
    def dfs(node, visited, length):
        nonlocal longest
        if length > longest:
            longest = length
        for nbr in graph[node]:
            if nbr not in visited:
                dfs(nbr, visited | {nbr}, length + 1)
                
    for idx in carbon_idxs:
        dfs(idx, {idx}, 1)
    
    return longest

def longest_chain_from_atom(mol, start, excluded):
    """
    Calculate the length of the longest chain (only counting carbon atoms) starting
    from the atom with index 'start' and not venturing into atoms in 'excluded'.
    """
    max_length = 0
    def dfs(current, visited, length):
        nonlocal max_length
        if length > max_length:
            max_length = length
        for nbr in current.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited and nbr.GetIdx() not in excluded:
                dfs(nbr, visited | {nbr.GetIdx()}, length + 1)
    dfs(mol.GetAtomWithIdx(start), {start}, 1)
    return max_length

def max_side_chain_length(mol, motif_atom_indices):
    """
    For a given motif (tuple of atom indices), look for alkyl side chains attached to any motif atom.
    Returns the maximum carbon chain length found among neighbors that are outside the motif.
    """
    max_chain = 0
    motif_set = set(motif_atom_indices)
    for idx in motif_atom_indices:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in motif_set:
                chain_len = longest_chain_from_atom(mol, nbr.GetIdx(), motif_set)
                if chain_len > max_chain:
                    max_chain = chain_len
    return max_chain

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid using improved heuristics.
    
    Heuristics:
      - Valid molecule with oxygen, at least 15 carbons and MW > 200 Da.
      - For phytocannabinoids: presence of a dihydroxybenzene (resorcinol-like) motif
        with an attached alkyl side chain of moderate length (5-10 carbons).
      - For synthetic cannabinoids: presence of an indole ring.
      - For endocannabinoids: a polar group (amide, ester or acid) plus an overall
        longest carbon chain of at least 16 atoms.
      - Fallback: an oxygenated heterocycle (ring of size ≥5 containing an oxygen)
        and a long overall carbon chain (≥16).
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if classified as cannabinoid, False if not.
      str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Preliminary checks.
    # Reject molecules containing phosphorus (not common in cannabinoids).
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Molecule contains phosphorus which is not typical in cannabinoid structures."
    
    # Must contain at least one oxygen.
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count == 0:
        return False, "Molecule lacks oxygen, but cannabinoids typically have oxygen in a ring or function."
    
    # Check total carbon count and molecular weight.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 15:
        return False, f"Carbon count ({total_carbons}) too low for a cannabinoid."
    
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight ({mol_wt:.1f} Da) too low for a cannabinoid."
    
    overall_chain = longest_carbon_chain(mol)

    # ----- Rule 1: Phytocannabinoids (Resorcinol-like motif with moderate alkyl side chain) -----
    # We use two SMARTS variants for resorcinol.
    resorcinol_patterns = [
        Chem.MolFromSmarts("c1cc(O)c(O)cc1"),
        Chem.MolFromSmarts("c1cc(O)cc(O)c1")
    ]
    for pattern in resorcinol_patterns:
        if pattern and mol.HasSubstructMatch(pattern):
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                side_chain = max_side_chain_length(mol, match)
                # Expect an alkyl side chain in the 5-10 carbon range (typical is about 8).
                if 5 <= side_chain <= 10:
                    return True, f"Molecule contains a resorcinol-like motif with an alkyl side chain (chain length {side_chain}), typical of phytocannabinoids."
    
    # ----- Rule 2: Synthetic cannabinoids (Indole ring) -----
    # Use an improved indole SMARTS.
    indole = Chem.MolFromSmarts("c1cc2c(c1)[nH]cc2")
    if indole and mol.HasSubstructMatch(indole):
        return True, "Molecule contains an indole ring, a characteristic feature of synthetic cannabinoids."
    
    # ----- Rule 3: Endocannabinoids (Polar group plus long carbon chain) -----
    # Define SMARTS for amide, ester, and carboxylic acid groups.
    amide  = Chem.MolFromSmarts("[NX3;!H0][CX3](=O)")
    ester  = Chem.MolFromSmarts("[CX3](=O)[OX2H0]")
    acid   = Chem.MolFromSmarts("[CX3](=O)[OX1H]")
    has_polar = False
    for pattern in (amide, ester, acid):
        if pattern and mol.HasSubstructMatch(pattern):
            has_polar = True
            break
    if has_polar and overall_chain >= 16:
        return True, f"Molecule has a polar group (amide/ester/acid) and a long carbon chain (chain length {overall_chain}), common features of endocannabinoids."
    
    # ----- Rule 4: Fallback: oxygenated heterocycle (ring size ≥5 with oxygen) with a long chain.
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) >= 5:
            # Check if the ring has an oxygen.
            if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring):
                if overall_chain >= 16:
                    return True, "Molecule contains an oxygenated heterocycle and a long carbon chain, possibly cannabinoid-related."
    
    return False, "Molecule does not exhibit key cannabinoid features (resorcinol/indole motif, or polar group with long chain, or suitable oxygenated heterocycle)."

# For testing purposes:
if __name__ == "__main__":
    test_smiles = [
        # Some examples from the list:
        "O=C(O[C@@H]([C@@H](O)[C@H](O)CO)CO)C(=CC(C(O)C(=CC(C(O)C(=CC(C(O[C@@H]1O[C@@H]([C@@H](O)[C@@H]([C@@H]1O)O)CO)C(=CC(CC(CC)C)C)C)C)C)C)C)C)C",  # Roselipin 3E; has an ester plus long chain → endocannabinoid
        "O(C(=O)CCCCCCCCC/C=C\\C/C=C\\CCCCC)C(CO)CO",  # MG(0:0/20:2(11Z,14Z)/0:0)
        "CCCC\\C=C/CCCCCCCC(=O)NCCO",  # palmitoleoyl ethanolamide
        "[H][C@]1(C=C(C)CC[C@H]1C1=C(O)C=C(CCCCC)C=C1O)C(C)=C",  # cannabidiol; resorcinol + side chain ~8
        "CCCCCN1C=C(C(=O)CC2=CC=C(OC)C=C2)C2=CC=CC=C12",  # JWH-201; synthetic cannabinoid (indole ring)
    ]
    
    for smi in test_smiles:
        result, reason = is_cannabinoid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")