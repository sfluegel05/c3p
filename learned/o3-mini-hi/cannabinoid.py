"""
Classifies: CHEBI:67194 cannabinoid
"""
"""
Classifies: Cannabinoid chemical entities

Definition:
  "A diverse group of pharmacologically active secondary metabolite characteristic to Cannabis plant as well as produced naturally in the body by humans and animals. 
   Cannabinoids contain oxygen as a part of the heterocyclic ring or in the form of various functional groups."

Heuristic rules in this version:
  1. Check basic requirements: valid molecule, contains oxygen, no unwanted elements (e.g. phosphorus), minimum carbon count (>=15)
     and molecular weight above ~200 Da.
  2. For phytocannabinoids: look for a dihydroxybenzene (resorcinol-like) motif using two different SMARTS patterns.
     Also require that at least one carbon outside the ring is connected to a long side chain (≥4 consecutive carbons).
  3. For synthetic cannabinoids: search for an indole ring.
  4. For endocannabinoids: if the molecule contains a polar group (amide, ester, or acid) AND the overall longest carbon chain is ≥16,
     it is classified.
  5. Fallback: if the molecule contains any oxygenated heterocyclic ring (of size ≥5) and the longest carbon chain is ≥16, call it cannabinoid‐related.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_carbon_chain(mol):
    """
    Returns the length (number of atoms) of the longest continuous chain of carbon atoms.
    We create a carbon–only connectivity graph and use DFS (Depth-first search) to search for the longest path.
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
    Calculate the length of the longest chain (counting only carbon atoms) starting
    from the atom with index 'start', without traversing any atom in 'excluded'.
    """
    max_length = 0
    def dfs(current, visited, length):
        nonlocal max_length
        if length > max_length:
            max_length = length
        for nbr in current.GetNeighbors():
            idx = nbr.GetIdx()
            if nbr.GetAtomicNum() == 6 and idx not in visited and idx not in excluded:
                dfs(nbr, visited | {idx}, length + 1)
    dfs(mol.GetAtomWithIdx(start), {start}, 1)
    return max_length

def max_side_chain_length(mol, motif_atom_indices):
    """
    For a given motif match (tuple of atom indices), examine each atom in the motif.
    For each neighbor (a carbon) not in the motif, compute the longest carbon chain length.
    Return the maximum side chain length found.
    """
    max_chain = 0
    motif_set = set(motif_atom_indices)
    for idx in motif_atom_indices:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            nid = nbr.GetIdx()
            if nbr.GetAtomicNum() == 6 and nid not in motif_set:
                chain_len = longest_chain_from_atom(mol, nid, motif_set)
                if chain_len > max_chain:
                    max_chain = chain_len
    return max_chain

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid using refined heuristics.
    
    Heuristics:
    - It must be a valid molecule with oxygen, sufficient carbon atoms (>=15) and adequate molecular weight (>200 Da).
    - For phytocannabinoids, look for a resorcinol-like motif (two SMARTS patterns) with at least one alkyl side chain (≥4 carbons) attached.
    - For synthetic cannabinoids, check for the presence of an indole ring.
    - For endocannabinoids, if a polar group (amide/ester/acid) is present along with an overall long carbon chain (≥16), classify as cannabinoid.
    - As a fallback, if the molecule contains an oxygenated heterocycle (≥5 atoms, at least one oxygen) and a long chain (≥16), classify it as potentially cannabinoid.
    
    Args:
        smiles (str): SMILES string representing the molecule.
        
    Returns:
        bool: True if the molecule is classified as a cannabinoid, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Exclude unwanted elements, e.g. phosphorus (atomic number 15).
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Molecule contains phosphorus which is not typical in cannabinoid structures"
    
    # Must contain oxygen.
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count == 0:
        return False, "Molecule does not contain oxygen, but cannabinoids typically do."
    
    # Check minimum number of carbons and molecular weight.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 15:
        return False, f"Total carbon count ({total_carbons}) is too low to be a typical cannabinoid."
    
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low to be a typical cannabinoid."
    
    # Compute the overall longest chain of carbons.
    overall_chain = longest_carbon_chain(mol)
    
    # ----- Rule 1: Phytocannabinoids (Resorcinol-like motif) -----
    # Use two SMARTS patterns to be flexible about the dihydroxybenzene arrangement.
    patterns = [
        Chem.MolFromSmarts("c1cc(O)c(O)cc1"),  # classic resorcinol
        Chem.MolFromSmarts("c1cc(O)cc(O)c1")   # alternative ordering
    ]
    for pattern in patterns:
        if pattern and mol.HasSubstructMatch(pattern):
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                side_chain = max_side_chain_length(mol, match)
                # Relax requirement slightly to a chain length of at least 4 carbons.
                if side_chain >= 4:
                    return True, f"Molecule contains a dihydroxybenzene (resorcinol-like) motif with an alkyl side chain (chain length {side_chain}), common in phytocannabinoids."
    
    # ----- Rule 2: Synthetic cannabinoids (Indole ring) -----
    indole = Chem.MolFromSmarts("c1ccc2[nH]ccc2c1")
    if indole and mol.HasSubstructMatch(indole):
        return True, "Molecule contains an indole ring, a feature common in synthetic cannabinoids."
    
    # ----- Rule 3: Endocannabinoids (Polar group plus long carbon chain) -----
    # Define SMARTS for amide, ester, and carboxylic acid groups.
    amide = Chem.MolFromSmarts("[NX3;!H0][CX3](=O)")
    ester = Chem.MolFromSmarts("[CX3](=O)[OX2H0]")
    acid  = Chem.MolFromSmarts("[CX3](=O)[OX1H]")
    has_polar = False
    for pattern in (amide, ester, acid):
        if pattern and mol.HasSubstructMatch(pattern):
            has_polar = True
            break
    if has_polar and overall_chain >= 16:
        return True, f"Molecule has a polar group (amide/ester/acid) and a long aliphatic chain (chain length {overall_chain}), features common in endocannabinoids."
    
    # ----- Rule 4: Fallback (Oxygenated heterocycle + long chain) -----
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) >= 5:
            if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring):
                if overall_chain >= 16:
                    return True, "Molecule contains an oxygenated heterocyclic ring and a long aliphatic chain, a possible cannabinoid feature."
    
    return False, "Molecule does not contain key cannabinoid structural features (resorcinol/indole core, nor polar group with long chain, nor oxygenated heterocycle with long chain)."

# For testing purposes:
if __name__ == "__main__":
    test_smiles = [
        "O=C(O[C@@H]([C@@H](O)[C@H](O)CO)CO)C(=CC(C(O)C(=CC(C(O)C(=CC(C(O[C@@H]1O[C@@H]([C@@H](O)[C@@H]([C@@H]1O)O)CO)C(=CC(CC(CC)C)C)C)C)C)C)C)C)C",  # Roselipin 3E
        "O(C(=O)CCCCCCCCC/C=C\\C/C=C\\CCCCC)C(CO)CO",  # MG(0:0/20:2(11Z,14Z)/0:0)
        "CCCC\\C=C/CCCCCCCC(=O)NCCO",  # palmitoleoyl ethanolamide
        "[H][C@]1(C=C(C)CC[C@H]1C1=C(O)C=C(CCCCC)C=C1O)C(C)=C",  # cannabidiol
        "CCCCCN1C=C(C(=O)CC2=CC=C(OC)C=C2)C2=CC=CC=C12",  # JWH-201
    ]
    
    for smi in test_smiles:
        result, reason = is_cannabinoid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")