"""
Classifies: CHEBI:27300 vitamin D
"""
"""
Classifies: vitamin D compounds (fat‐soluble hydroxy seco‐steroids)

Heuristic criteria (revised):
  1. SMILES must be valid.
  2. Molecular weight roughly in the range 250–700 Da.
  3. Must contain at least one hydroxyl group (-OH).
  4. Must contain a “seco‐steroid” core. For our purposes:
       – Consider only rings of size 5 or 6 that are mainly carbocyclic (≥4 carbons).
       – Build a graph over these rings; then allow the largest fused component to be either:
             • exactly 3 rings with sorted sizes [5,6,6], or
             • exactly 2 rings with sorted sizes [5,6] or [6,6] (to allow for the “seco” opening).
  5. There must be at least one exocyclic double bond (a C=C where at least one atom is not in a ring).
  6. The molecule should be lipophilic (calculated logP ≥ 3.0).

If any test fails, the function returns False plus a reason.
Note: This heuristic filter may misclassify borderline cases.
"""
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D compound (a fat‐soluble hydroxy seco‐steroid)
    based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule meets vitamin D heuristic criteria, False otherwise.
        str: Explanation for the decision.
    """
    # 1. Parse SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 2. Molecular weight check: roughly 250–700 Da.
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 700:
        return False, f"Molecular weight {mol_wt:.1f} Da is outside the typical range (250–700 Da) for vitamin D compounds"
    
    # 3. Must contain at least one hydroxyl group (-OH).
    hydroxyl_smarts = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(hydroxyl_smarts):
        return False, "No hydroxyl (-OH) group found; vitamin D compounds must have at least one -OH"
    
    # 4. Identify candidate rings for the steroid core.
    # Consider rings of size 5 or 6 that are largely carbocyclic (≥4 carbons).
    ri = mol.GetRingInfo()
    all_rings = ri.AtomRings()  # each ring is a tuple of atom indices
    candidate_rings = []
    ring_sizes = []  # will store size for each candidate ring
    for ring in all_rings:
        if len(ring) in (5, 6):
            n_carbons = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            if n_carbons >= 4:
                candidate_rings.append(set(ring))
                ring_sizes.append(len(ring))
    
    if not candidate_rings:
        return False, "No 5- or 6-membered (mostly carbocyclic) rings found; lacking a steroid-like core"
    
    # Build an adjacency graph of candidate rings: two rings are connected if they share at least one atom.
    n = len(candidate_rings)
    adj = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            if candidate_rings[i].intersection(candidate_rings[j]):
                adj[i].add(j)
                adj[j].add(i)
    
    # Use depth-first search (DFS) to find connected components (fused ring systems)
    visited = set()
    components = []
    def dfs(node, comp):
        comp.add(node)
        for neigh in adj[node]:
            if neigh not in comp:
                dfs(neigh, comp)
    for i in range(n):
        if i not in visited:
            comp = set()
            dfs(i, comp)
            visited |= comp
            components.append(comp)
    
    if not components:
        return False, "No fused ring system found; lacking a seco-steroid core"
    
    # Examine the largest fused set.
    largest_comp = max(components, key=lambda comp: len(comp))
    n_rings = len(largest_comp)
    # Allow either a 3-ring fused set or (if the B-ring is cleaved) a 2-ring system.
    if n_rings not in (2, 3):
        return False, f"Fused ring system has {n_rings} rings; expected 2 (seco) or 3 rings for a vitamin D core"
    
    # Get the ring sizes (sorted) for rings in the largest fused component.
    comp_ring_sizes = sorted([len(list(candidate_rings[i])) for i in largest_comp])
    if n_rings == 3:
        if comp_ring_sizes != [5, 6, 6]:
            return False, f"Fused 3-ring system sizes are {comp_ring_sizes}; expected [5, 6, 6] for a vitamin D core"
    elif n_rings == 2:
        # Accept [6,6] or [5,6] patterns when one ring is cleaved.
        if comp_ring_sizes not in ([6, 6], [5, 6]):
            return False, f"Fused 2-ring system sizes are {comp_ring_sizes}; expected [6,6] or [5,6] for a vitamin D core (seco-steroid)"
    
    # 5. Check for at least one exocyclic double bond.
    # A bond is exocyclic if it is a double bond and at least one atom is not in any ring.
    exocyclic_db_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if not (a1.IsInRing() and a2.IsInRing()):
                exocyclic_db_count += 1
    if exocyclic_db_count < 1:
        return False, "No exocyclic double bond found; vitamin D compounds typically have an unsaturated side chain"
    
    # 6. Lipophilicity (logP) check.
    logP = Crippen.MolLogP(mol)
    if logP < 3.0:
        return False, f"Calculated logP ({logP:.1f}) is too low; vitamin D compounds are expected to be lipophilic (logP ≥ 3.0)"
    
    return True, ("Molecule meets vitamin D criteria: appropriate molecular weight, contains at least one -OH, "
                   "a seco-steroid core (fused rings with acceptable sizes), an exocyclic double bond, "
                   "and sufficient lipophilicity.")

# For testing when running as a standalone script.
if __name__ == "__main__":
    # Example vitamin D SMILES from the provided examples:
    test_smiles_list = [
        "C[C@H](C=C[C@H](C)C(C)C)C1CC[C@@H]2[C@@]1(CCCC2=CC=C3C[C@H](CCC3=C)O)C",  # (1S)-3-[...]-4-cyclohexanol
        "O([C@@H](C[C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\\3/C[C@@H](O)C[C@H](O)C3=C)[H])C)[H])C)C(=O)C(O)(C)C)C4O[C@@H]([C@@H](O)[C@H](O)[C@H]4O)C(O)=O"  # (23S)-23,25-dihdroxy-... glucuronide
    ]
    
    for smi in test_smiles_list:
        result, reason = is_vitamin_D(smi)
        print(result, reason)