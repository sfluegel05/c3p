"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
"""
Classifies: CHEBI:16β–hydroxy steroid
Definition: A 16–hydroxy steroid in which the hydroxy group at position 16 has a beta–configuration.
Heuristic improvements over the previous attempt:
  1. We again extract rings of size 5 or 6 and group fused rings.
  2. For each connected component with ≥4 rings we examine all sets of exactly 4 fused rings.
  3. We now require that among these four rings there is exactly one 5–membered ring and exactly three 6–membered rings.
  4. Critically, we require that the union of atoms in the four rings (i.e. the steroid nucleus) contains exactly 17 carbon atoms.
  5. Then—in the candidate D–ring (the only 5–membered ring)—we require that at least one carbon atom, which has defined stereochemistry,
     is directly bound (by a single bond) to an oxygen atom that appears to be a hydroxyl group (i.e. with at least one hydrogen).
     
NOTE: This heuristic still does not capture all steroid variants perfectly but it is tuned to reduce false positives.
"""
from rdkit import Chem
from rdkit.Chem import rdchem
from itertools import combinations

def is_16beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 16β–hydroxy steroid via an improved heuristic.
    The molecule must possess a fused four–ring system (steroid nucleus) composed of three six–membered rings 
    and one five–membered ring. The union of atoms in these rings must yield exactly 17 carbon atoms.
    In the unique five–membered ring (candidate D–ring) at least one carbon with defined stereochemistry 
    must have a single bond to an oxygen that appears to be an –OH group (has at least one hydrogen).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule qualifies as a 16β–hydroxy steroid (by our heuristic), False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First ensure stereochemistry is perceived
    Chem.AssignStereochemistry(mol, force=True)
    
    # Retrieve all ring systems from the molecule and restrict to rings of size 5 or 6.
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()  # Each is a tuple of atom indices.
    rings_5_6 = []
    for ring in all_rings:
        if len(ring) in (5, 6):
            rings_5_6.append(set(ring))
    if not rings_5_6:
        return False, "No rings of size 5 or 6 found; cannot be a steroid nucleus"
    
    # Build a graph where nodes represent rings (from rings_5_6) and edges connect rings sharing >=2 atoms.
    num_rings = len(rings_5_6)
    ring_graph = {i: set() for i in range(num_rings)}
    for i in range(num_rings):
        for j in range(i+1, num_rings):
            if len(rings_5_6[i].intersection(rings_5_6[j])) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Find connected components of rings.
    seen = set()
    components = []
    for i in range(num_rings):
        if i not in seen:
            stack = [i]
            comp = set()
            while stack:
                current = stack.pop()
                if current not in seen:
                    seen.add(current)
                    comp.add(current)
                    stack.extend(ring_graph[current])
            components.append(comp)
    
    candidate_nucleus = None
    candidate_D_ring = None  # this is the unique 5–membered ring in the nucleus.
    
    # Look in each connected component for a combination of exactly 4 rings with the required properties.
    for comp in components:
        if len(comp) < 4:
            continue
        comp_list = list(comp)
        for comb in combinations(comp_list, 4):
            rings_subset = [rings_5_6[i] for i in comb]
            # Check ring size distribution: exactly 1 five-membered and 3 six-membered rings.
            count5 = sum(1 for ring in rings_subset if len(ring) == 5)
            count6 = sum(1 for ring in rings_subset if len(ring) == 6)
            if count5 != 1 or count6 != 3:
                continue
            # Verify connectivity: each ring in the subset must be fused (share >=2 atoms) with at least one other ring.
            mini_connected = True
            for i in range(4):
                if not any(i != j and len(rings_subset[i].intersection(rings_subset[j])) >= 2 
                           for j in range(4)):
                    mini_connected = False
                    break
            if not mini_connected:
                continue
            # Get the union of all atoms in the candidate nucleus.
            union_atoms = set()
            for ring in rings_subset:
                union_atoms.update(ring)
            # Count carbon atoms among these atoms.
            c_count = sum(1 for idx in union_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            # For a classical steroid nucleus (cyclopentanoperhydrophenanthrene) there should be exactly 17 carbons.
            if c_count != 17:
                continue
            
            # Identify the candidate D–ring (the only 5-membered ring).
            for ring in rings_subset:
                if len(ring) == 5:
                    # Require that the 5-membered ring is fused with at least one six–membered ring (by >=2 atoms).
                    if not any((len(ring.intersection(other)) >= 2) for other in rings_subset if other is not ring and len(other) == 6):
                        continue
                    candidate_D_ring = ring
                    candidate_nucleus = union_atoms
                    break
            if candidate_nucleus is not None and candidate_D_ring is not None:
                break
        if candidate_nucleus is not None:
            break
    
    if candidate_nucleus is None:
        return False, "No fused 4–ring steroid nucleus (three 6–membered and one 5–membered rings with exactly 17 carbons) found"
    
    # Within the candidate D–ring, look for a carbon (with defined chiral tag) that is directly bonded (single bond)
    # to an oxygen atom that appears to be in an –OH group (has at least one hydrogen).
    found_16beta_oh = False
    for idx in candidate_D_ring:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() != 6:
            continue
        # Exclude atoms without a specified chiral center.
        if atom.GetChiralTag() == rdchem.ChiralType.CHI_UNSPECIFIED:
            continue
        # Check neighbors for an oxygen bound by a single bond.
        for nb in atom.GetNeighbors():
            if nb.GetAtomicNum() != 8:
                continue
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx())
            if bond is None or bond.GetBondType() != rdchem.BondType.SINGLE:
                continue
            if nb.GetTotalNumHs() < 1:
                continue
            # If we reach here, this atom represents a chiral center in the D–ring with an –OH.
            found_16beta_oh = True
            break
        if found_16beta_oh:
            break
    
    if not found_16beta_oh:
        return False, "Steroid nucleus found but no chiral carbon with a directly attached –OH in the D–ring was detected"
    
    return True, ("Molecule contains a fused steroid nucleus (three 6–membered and one 5–membered rings with exactly 17 carbons) "
                  "and the candidate D–ring has a chiral carbon bearing an –OH group (heuristic for 16β–hydroxy steroid)")

# Example usage:
if __name__ == '__main__':
    # Test example: 16beta-hydroxytestosterone
    test_smiles = "C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1C[C@H](O)[C@@H]2O"
    result, reason = is_16beta_hydroxy_steroid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)