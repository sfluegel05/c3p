"""
Classifies: CHEBI:27300 vitamin D
"""
"""
Classifies: vitamin D compounds (fat‐soluble hydroxy seco‐steroids).

Heuristic criteria (revised):
  1. SMILES must be valid.
  2. Molecular weight roughly in the range 250–700 Da.
  3. Must contain at least one hydroxyl group (-OH).
  4. Must contain a "seco‐steroid" core. For our purposes we:
       – Consider only rings of size 5 or 6 that are largely carbocyclic (≥4 carbons).
       – Build a graph of rings that share atoms and require that the largest connected
         (fused) component comprises exactly 3 rings.
       – In addition, we check that the sorted ring sizes are exactly [5,6,6].
  5. Instead of a full conjugated diene motif we now require that the molecule has at 
     least one exocyclic double bond (i.e. a C=C bond for which at least one of the atoms is not in any ring).
  6. Shows significant lipophilicity (calculated logP ≥ 3.0).
If any test fails, the function returns a corresponding reason.
Note: This heuristic filter may misclassify borderline cases.
"""
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D compound (fat‐soluble hydroxy seco‐steroid)
    based on its SMILES string.
    
    Heuristic criteria:
      - Valid SMILES.
      - Molecular weight roughly within 250 to 700 Da.
      - Contains at least one hydroxyl group (-OH).
      - Contains a secosteroid core defined as a fused ring system among 5- or 6-membered rings
        (with at least 4 carbons per ring) where the largest connected group has exactly 3 rings
        and where the collection of ring sizes (sorted) is [5, 6, 6].
      - Contains at least one exocyclic double bond (a C=C where at least one atom is not in any ring).
      - Has sufficient lipophilicity (calculated logP ≥ 3.0).
      
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule meets vitamin D heuristic criteria, False otherwise.
        str: Explanation for the classification decision.
    """
    # 1. Validate and create RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 2. Molecular weight check (approx. 250–700 Da).
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 700:
        return False, f"Molecular weight {mol_wt:.1f} Da is outside the typical range (250–700 Da) for vitamin D compounds"
    
    # 3. Check for at least one hydroxyl (-OH) group.
    hydroxyl_smarts = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(hydroxyl_smarts):
        return False, "No hydroxyl (-OH) group found; vitamin D compounds must have at least one -OH"
    
    # 4. Identify candidate rings for the secosteroid core.
    # We consider only rings of size 5 or 6 that are largely carbocyclic (≥4 carbon atoms in the ring).
    ri = mol.GetRingInfo()
    all_rings = ri.AtomRings()  # Each ring is a tuple of atom indices.
    core_rings = []
    ring_sizes = []  # Save the size for later evaluation.
    for ring in all_rings:
        if len(ring) in (5, 6):
            n_carbons = sum(1 for i in ring if mol.GetAtomWithIdx(i).GetAtomicNum() == 6)
            if n_carbons >= 4:
                core_rings.append(set(ring))
                ring_sizes.append(len(ring))
    if not core_rings:
        return False, "No 5- or 6-membered (mostly carbocyclic) rings found; lacking a steroid-like core"
    
    # Build an adjacency graph over our candidate rings: each ring node shares an edge if the rings share at least one atom.
    n = len(core_rings)
    adj = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            if core_rings[i].intersection(core_rings[j]):
                adj[i].add(j)
                adj[j].add(i)
    
    # Find connected (fused) ring systems using DFS.
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
    
    # If no connected set found then reject.
    if not components:
        return False, "No fused ring system found; lacking a secosteroid core"
    
    # We require that the largest fused system has exactly 3 rings.
    largest_comp = max(components, key=lambda comp: len(comp))
    if len(largest_comp) != 3:
        return False, f"Fused ring system has {len(largest_comp)} rings; expected exactly 3 rings (a 5/6/6 pattern) for a vitamin D (seco-steroid) core"
    
    # Obtain the ring sizes of the rings in the largest fused component.
    comp_ring_sizes = sorted([len(list(core_rings[i])) for i in largest_comp])
    if comp_ring_sizes != [5, 6, 6]:
        return False, f"Fused ring system sizes are {comp_ring_sizes}; expected a 5,6,6 pattern for a vitamin D secosteroid core"
    
    # 5. Check for exocyclic double bonds.
    # We count a bond as exocyclic if it is a double bond and at least one atom is not in any ring.
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
        return False, f"Calculated logP ({logP:.1f}) is rather low; vitamin D compounds are expected to be lipophilic (logP ≥ 3.0)"
    
    return True, ("Molecule meets vitamin D criteria: appropriate molecular weight, contains at least one -OH, "
                   "a secosteroid core with a 5,6,6 fused ring system, an unsaturated (exocyclic double bond) side chain, "
                   "and sufficient lipophilicity.")

# Example usage when running this module directly.
if __name__ == "__main__":
    # Test with one vitamin D compound from the examples.
    test_smiles = "C[C@H](C=C[C@H](C)C(C)C)C1CC[C@@H]2[C@@]1(CCCC2=CC=C3C[C@H](CCC3=C)O)C"
    result, reason = is_vitamin_D(test_smiles)
    print(result, reason)