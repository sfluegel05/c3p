"""
Classifies: CHEBI:16158 steroid sulfate
"""
"""
Classifies: steroid sulfate
Definition: A sulfuric ester obtained by the formal condensation of a hydroxy group 
of any steroid with sulfuric acid. For our purposes a valid steroid sulfate must have a sulfate ester group (either anionic or neutral) 
that is attached either directly or via one intervening saturated (sp3 non‐aromatic) linker to a steroid nucleus.
A steroid nucleus is heuristically defined as a fused polycyclic system composed solely of 5‐ and 6‐membered rings,
with at least four rings overall (and with at least one 5‐membered ring and three 6‐membered rings). 
The improvement over the previous attempt is that now “fused” rings are only merged when they share an entire bond (≥2 atoms)
and the connectivity from the sulfate oxygen to the nucleus is checked via a shortest‐path search.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_steroid_sulfate(smiles: str):
    """
    Determines if a molecule is a steroid sulfate based on its SMILES string.
    
    Process:
    1. Parse the molecule.
    2. Identify candidate steroid nucleus:
       - Collect rings of size 5 or 6.
       - Build a graph of rings that are truly fused (share ≥2 atoms).
       - For each connected group, count the rings and their sizes.
         Accept if one group contains at least 4 rings (≥1 five‐membered and ≥3 six‐membered).
    3. Identify candidate sulfate group(s):
       - For each oxygen atom bonded to a sulfur via a single bond,
         check that the connected sulfur is “sulfate‐like” (bonded to at least 3 oxygens via single bonds).
    4. For each candidate sulfate ester oxygen, look for a connection to the nucleus:
       - Compute the shortest bond path from the candidate oxygen to any nucleus atom.
         Accept if the path length is 1 (direct attachment).
       - Or if the path length is 2, check that the intervening atom is sp3 and non‐aromatic.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a steroid sulfate; False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # ----- Step 1: Identify candidate steroid nucleus -----
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    if not all_rings:
        return False, "No rings found; not a steroid nucleus"
    
    # Consider only rings of size 5 or 6.
    rings_5_or_6 = [set(ring) for ring in all_rings if len(ring) in (5, 6)]
    if not rings_5_or_6:
        return False, "No 5- or 6-membered rings found; not a steroid nucleus"
    
    # Build a graph among rings - nodes are indices of rings in rings_5_or_6;
    # add an edge if two rings share at least 2 atoms (an entire bond).
    n_rings = len(rings_5_or_6)
    adj = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        for j in range(i+1, n_rings):
            # If the two rings share at least 2 atoms, consider them fused.
            if len(rings_5_or_6[i] & rings_5_or_6[j]) >= 2:
                adj[i].add(j)
                adj[j].add(i)
    
    # Determine connected components in the ring graph.
    visited = set()
    fused_components = []
    for i in range(n_rings):
        if i in visited:
            continue
        comp = set()
        stack = [i]
        while stack:
            node = stack.pop()
            if node in visited:
                continue
            visited.add(node)
            comp.add(node)
            stack.extend(adj[node] - visited)
        fused_components.append(comp)
    
    # For each fused component, get the union of atoms, and count how many rings (by original list)
    steroid_nucleus = None
    for comp in fused_components:
        # comp is a set of indices into rings_5_or_6
        component_rings = [rings_5_or_6[i] for i in comp]
        # Count rings by size (note that rings may share atoms, but count the ring instances)
        count_5 = sum(1 for ring in component_rings if len(ring) == 5)
        count_6 = sum(1 for ring in component_rings if len(ring) == 6)
        if len(component_rings) >= 4 and count_5 >= 1 and count_6 >= 3:
            # Use the union of all atoms in these rings as the nucleus
            nucleus_atoms = set()
            for ring in component_rings:
                nucleus_atoms.update(ring)
            steroid_nucleus = nucleus_atoms
            break
    
    if steroid_nucleus is None:
        return False, ("No fused tetracyclic (steroid) nucleus found "
                       "(requires at least 4 fused rings [with ≥1 five-membered and ≥3 six-membered rings])")
    
    # ----- Step 2: Identify candidate sulfate ester oxygen(s) -----
    candidate_oxy_indices = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue  # only oxygens
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 16:  # sulfur
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                # require single bond (to exclude S=O double bonds)
                if bond.GetBondType() != rdchem.BondType.SINGLE:
                    continue
                # Check that this sulfur is “sulfate‐like”: at least three oxygen neighbors via single bonds.
                oxy_count = 0
                for snbr in nbr.GetNeighbors():
                    if snbr.GetAtomicNum() == 8:
                        bnd = mol.GetBondBetweenAtoms(nbr.GetIdx(), snbr.GetIdx())
                        if bnd.GetBondType() == rdchem.BondType.SINGLE:
                            oxy_count += 1
                if oxy_count < 3:
                    continue
                # We have a candidate connection oxygen for the sulfate group.
                candidate_oxy_indices.append(atom.GetIdx())
                break  # no need to check further neighbors for this oxygen

    if not candidate_oxy_indices:
        return False, "No sulfate ester attachment oxygen found"

    # ----- Step 3: Check connectivity from candidate sulfate oxygen to the steroid nucleus -----
    # For each candidate oxygen, we check the shortest path (number of bonds) to ANY nucleus atom.
    attached_found = False
    for oxy_idx in candidate_oxy_indices:
        for nuc_idx in steroid_nucleus:
            # Get shortest path between candidate oxygen and this nucleus atom.
            path = Chem.GetShortestPath(mol, oxy_idx, nuc_idx)
            if not path:
                continue
            path_len = len(path) - 1  # number of bonds
            # Direct attachment:
            if path_len == 1:
                attached_found = True
                break
            # One-linker attachment is allowed if path length is 2.
            if path_len == 2:
                # path structure: [candidate oxygen] - X - nucleus_atom.
                linker_idx = path[1]
                linker = mol.GetAtomWithIdx(linker_idx)
                if linker.GetHybridization() == rdchem.HybridizationType.SP3 and not linker.GetIsAromatic():
                    attached_found = True
                    break
        if attached_found:
            break

    if not attached_found:
        return False, "Sulfate group not attached (even via one saturated sp3 linker) to the steroid nucleus"

    return True, ("Contains sulfate ester group attached (directly or via one sp³ linker) "
                  "to steroid nucleus (fused tetracyclic system)")

# Example usage (testing with sample SMILES):
if __name__ == "__main__":
    test_smiles = [
        # True positives:
        "[Na+].C[C@]12CC[C@H]3C(=CCc4cc(OS([O-])(=O)=O)ccc34)[C@@H]1CCC2=O",  # equilin sodium sulfate
        "C[C@]12CC[C@@H]3[C@H]([C@H]1CCC2=O)CCC4=C3C=CC(=C4)OS(=O)(=O)O",  # sulfuric acid ester example
        "S(O[C@H]1C[C@@]2([C@@]3([C@]([C@]4([C@](CC3)([C@H]([C@@H](O)[C@H]4O)C=C)C)[H])(CC[C@]2([C@@H](O)[C@H]1O)[H])[H])[H])C)(O)(=O)=O",  # Ptilosteroid A
        "S(O[C@H]1C[C@@]2([C@@](C3C(C4[C@](CC3)(C(=O)CC4)C)CC2)(CC1)C)[H])(O)(=O)=O",  # Etiocholanolone sulfate
        # False negatives (should be classified as steroid sulfate):
        "S(OCC(=O)[C@@H]1[C@@]2([C@]([C@]3([C@@]([C@@]4(C(=CC3)C[C@@H](O)CC4)C)(CC2)[H])[H])(CC1)[H])C)(O)(=O)=O",  # 21-hydroxypregnenolone monosulfate
        # One more: the examples (e.g., 5beta-scymnol sulfate) where the sulfate is connected via a linker.
        "[H][C@@]12C[C@H](O)CC[C@]1(C)[C@@]1([H])C[C@H](O)[C@]3(C)[C@]([H])(CC[C@@]3([H])[C@]1([H])[C@H](O)C2)[C@H](C)CCCC(CO)COS(O)(=O)=O",
        # Molecules that are NOT steroid sulfates:
        "O[S](=O)(=O)O",  # pure sulfate group
        "CCCC",         # non-cyclic aliphatic
    ]
    for smi in test_smiles:
        result, reason = is_steroid_sulfate(smi)
        print("SMILES:", smi)
        print("Result:", result, "|", reason)
        print("-"*80)