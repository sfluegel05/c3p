"""
Classifies: CHEBI:67194 cannabinoid
"""
"""
Classifies: Cannabinoid chemical entities

Definition:
  "A diverse group of pharmacologically active secondary metabolite characteristic to Cannabis plant as well as produced naturally in the body by humans and animals. 
   Cannabinoids contain oxygen as a part of the heterocyclic ring or in the form of various functional groups."

Heuristic rules in this improved version:
  1. Molecule must be valid, contain oxygen, have no unwanted elements (e.g. phosphorus), have at least 15 carbons and a molecular weight above ~200 Da.
  2. For phytocannabinoids: look for a dihydroxybenzene (resorcinol) motif (SMARTS "c1cc(O)c(O)cc1") and require that at least one carbon outside the ring is attached as a long side chain (≥5 consecutive carbons).
  3. For synthetic cannabinoids: check for an indole ring using an alternate SMARTS "c1ccc2[nH]ccc2c1".
  4. For endocannabinoids: if the molecule contains a polar group (ester, amide, or carboxylic acid) AND the overall longest carbon chain is ≥16, it is classified.
  5. As a fallback, if the molecule contains a heterocyclic ring (≥5 atoms) that includes oxygen and the molecule has an overall long carbon chain (≥16) then it might be cannabinoid‐related.
  
If none of these rules apply then the molecule is not classified as a cannabinoid.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_carbon_chain(mol):
    """
    Returns the length (number of atoms) of the longest chain
    of carbon atoms in the molecule. We build a simple carbon–only connectivity graph
    and use DFS to search for the longest simple path.
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
    Calculate the length of the longest chain (in number of carbon atoms) starting
    from a given atom index, without traversing any atom in 'excluded'.
    Only carbons (atomic number 6) are considered.
    """
    max_length = 0
    def dfs(current, visited, length):
        nonlocal max_length
        # update maximum chain length found so far
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
    For a given match (tuple of atom indices for the motif), check all atoms in the motif.
    For each neighbor (which is a carbon) not in the motif, compute the longest chain length.
    Return the maximum chain length found.
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

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a cannabinoid, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Exclude unwanted elements such as phosphorus (atomic number 15).
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Molecule contains phosphorus which is not typical in cannabinoid structures"
    
    # Must contain oxygen.
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count == 0:
        return False, "Molecule does not contain oxygen, but cannabinoids typically require oxygen."
    
    # Minimum carbon count and molecular weight.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 15:
        return False, f"Total carbon count ({total_carbons}) is too low to be a typical cannabinoid."
    
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low to be a typical cannabinoid."
    
    # Overall longest carbon chain.
    overall_chain = longest_carbon_chain(mol)
    
    # Rule 1 (Phytocannabinoids): Check for dihydroxybenzene (resorcinol) motif.
    resorcinol = Chem.MolFromSmarts("c1cc(O)c(O)cc1")
    if resorcinol and mol.HasSubstructMatch(resorcinol):
        matches = mol.GetSubstructMatches(resorcinol)
        for match in matches:
            side_chain = max_side_chain_length(mol, match)
            # Require that at least one alkyl chain of length >=5 is attached to the dihydroxybenzene ring.
            if side_chain >= 5:
                return True, f"Molecule contains a resorcinol (dihydroxybenzene) motif with an alkyl side chain (chain length {side_chain}), common in phytocannabinoids."
        # If motif found but no long side chain, do not classify as cannabinoid based solely on this rule.
    
    # Rule 2 (Synthetic cannabinoids): Look for an indole ring.
    indole = Chem.MolFromSmarts("c1ccc2[nH]ccc2c1")
    if indole and mol.HasSubstructMatch(indole):
        return True, "Molecule contains an indole ring, a feature common in synthetic cannabinoids."
    
    # Rule 3 (Endocannabinoids): Check for polar groups paired with a long aliphatic chain.
    # Define SMARTS for amide, ester, and carboxylic acid groups.
    amide = Chem.MolFromSmarts("[NX3!H0][CX3](=O)")
    ester = Chem.MolFromSmarts("[CX3](=O)[OX2H0]")
    acid  = Chem.MolFromSmarts("[CX3](=O)[OX1H]")
    has_polar = any(mol.HasSubstructMatch(pattern) for pattern in (amide, ester, acid) if pattern is not None)
    if has_polar and overall_chain >= 16:
        return True, f"Molecule has a polar group (amide/ester/acid) and a long aliphatic chain (chain length {overall_chain}), features common in endocannabinoids."
    
    # Rule 4 (Fallback): Check for any heterocyclic ring (size ≥5) that includes at least one oxygen,
    # combined with an overall long carbon chain.
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) >= 5:
            if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring):
                if overall_chain >= 16:
                    return True, "Molecule contains a heterocyclic ring with oxygen and a long aliphatic chain, a possible cannabinoid feature."
    
    # If none of the rules apply, do not classify.
    return False, "Molecule does not contain key cannabinoid structural features (resorcinol/indole core, or polar group with long chain, or oxygenated heterocycle with long chain)."

# Example calls (for testing; uncomment to run):
# test_smiles = [
#    "O=C(O[C@@H]([C@@H](O)[C@H](O)CO)CO)C(=CC(C(O)C(=CC(C(O)C(=CC(C(O[C@@H]1O[C@@H]([C@@H](O)[C@@H]([C@@H]1O)O)CO)C(=CC(CC(CC)C)C)C)C)C)C)C)C)C",  # Roselipin 3E
#    "O(C(=O)CCCCCCCCC/C=C\\C/C=C\\CCCCC)C(CO)CO",  # MG(0:0/20:2(11Z,14Z)/0:0)
#    "CCCC\\C=C/CCCCCCCC(=O)NCCO",  # palmitoleoyl ethanolamide
#    "[H][C@]1(CCC(C)=C[C@H]1C1=C(O)C=C(CCCCC)C=C1O)C(C)=C",  # cannabidiol
#    "CCCCCN1C=C(C(=O)CC2=CC=C(OC)C=C2)C2=CC=CC=C12",  # JWH-201
# ]
#
# for smi in test_smiles:
#     result, reason = is_cannabinoid(smi)
#     print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")