"""
Classifies: CHEBI:67194 cannabinoid
"""
"""
Classifies: Cannabinoid class
Definition:
  "A diverse group of pharmacologically active secondary metabolite characteristic to Cannabis plant 
   as well as produced naturally in the body by humans and animals. Cannabinoids contain oxygen as a 
   part of the heterocyclic ring or in the form of various functional groups."
   
Heuristic rules used in this improved version:
  1. Molecule must be a valid structure with at least one oxygen and a molecular weight above a minimum (e.g., 200 Da).
  2. Molecules with phosphorus (or similar non-cannabinoid markers) are excluded.
  3. If the molecule contains an aromatic ring with a hydroxyl group (phenol), then it is likely a phytocannabinoid.
  4. Alternatively, if the molecule contains an amide or ester functionality and features a sufficiently long aliphatic (carbon) chain 
     (using a DFS search to compute longest chain; here a chain of at least 8 carbons is required), it is classified as an endocannabinoid.
  5. As a last resort, if the molecule contains a heterocyclic ring (ring of at least 5 atoms) containing oxygen, then it might be cannabinoid-related.
     
If none of these rules apply, then the molecule is not classified as a cannabinoid.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_carbon_chain(mol):
    """
    Returns the length (number of atoms) of the longest chain 
    of carbon atoms in the molecule. It builds an undirected graph
    of carbon atoms (atomic number==6) and finds the longest simple path.
    """
    # Get all carbon atom indices
    carbon_idxs = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    # Build a connectivity graph among carbon atoms only.
    graph = {idx: [] for idx in carbon_idxs}
    for bond in mol.GetBonds():
        a = bond.GetBeginAtomIdx()
        b = bond.GetEndAtomIdx()
        if a in graph and b in graph:
            graph[a].append(b)
            graph[b].append(a)
    
    longest = 0
    # Use recursion to explore all simple paths starting from a given carbon
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

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on selected heuristics.
    
    Approach:
      - Validate the molecule and ensure it contains oxygen and avoids elements like phosphorus.
      - If the molecule is too small (low molecular weight), rule it out.
      - If an aromatic ring with a hydroxyl group (phenol) is found, assume the molecule is a phytocannabinoid.
      - Alternatively, if an amide or ester functionality exists together with a long aliphatic chain (longest C–chain length >= 8),
        then assume it is an endocannabinoid.
      - Otherwise, if any heterocyclic ring (ring of >=5 atoms that contains oxygen) is found, then return a positive classification.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a cannabinoid, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Exclude molecules with phosphorus (P atomic number 15), which are unlikely cannabinoids.
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Molecule contains phosphorus which is not typical in cannabinoid structures"
    
    # Check that the molecule actually contains oxygen.
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count == 0:
        return False, "Molecule does not contain any oxygen atoms, but cannabinoids require oxygen."
    
    # Check molecular weight (using exact weight). Cannabinoids commonly have MW above ~200 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low to be a typical cannabinoid."
    
    # Rule 1: Look for a phenol substructure:
    # Using a SMARTS pattern for an aromatic ring with a hydroxyl group.
    phenol_pattern = Chem.MolFromSmarts("c1ccc(O)cc1")
    if phenol_pattern and mol.HasSubstructMatch(phenol_pattern):
        return True, "Molecule contains an aromatic ring with a hydroxyl group (phenol), a common feature of phytocannabinoids."
    
    # Rule 2: Look for polar groups (amide or ester) using SMARTS patterns.
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)")
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    has_amide = mol.HasSubstructMatch(amide_pattern) if amide_pattern else False
    has_ester = mol.HasSubstructMatch(ester_pattern) if ester_pattern else False
    
    # Use DFS based longest carbon chain calculation (heuristic for long aliphatic chain).
    chain_length = longest_carbon_chain(mol)
    
    if (has_amide or has_ester) and chain_length >= 8:
        return True, f"Molecule has a polar head group (amide/ester) and a long aliphatic chain (chain length {chain_length}), features common in endocannabinoids."
    
    # Rule 3: Check if the molecule contains a heterocyclic ring with oxygen.
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) >= 5:
            # If any atom in the ring is oxygen, we consider it a heterocycle of interest.
            if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring):
                return True, "Molecule contains a heterocyclic ring with oxygen, a possible cannabinoid feature."
    
    return False, "Molecule does not contain key cannabinoid structural features (phenolic aromatic ring, or polar head group with long chain, or heterocyclic ring with oxygen)."

# Example calls (uncomment to test):
# smiles_examples = [
#     "O=C(O[C@@H]([C@@H](O)[C@H](O)CO)CO)C(=CC(C(O)C(=CC(C(O)C(=CC(C(O[C@@H]1O[C@@H]([C@@H](O)[C@@H]([C@@H]1O)O)CO)C(=CC(CC(CC)C)C)C)C)C)C)C)C)C", # Roselipin 3E
#     "O(C(=O)CCCCCCCCC/C=C\\C/C=C\\CCCCC)C(CO)CO", # MG(0:0/20:2(11Z,14Z)/0:0)
#     "CCCC\C=C/CCCCCCCC(=O)NCCO", # palmitoleoyl ethanolamide
# ]
# for smi in smiles_examples:
#     result, reason = is_cannabinoid(smi)
#     print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")