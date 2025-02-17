"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
"""
Classifies: CHEBI:16β–hydroxy steroid
Definition: A 16–hydroxy steroid in which the hydroxy group at position 16 has a beta–configuration.
Heuristic improvements:
  1. Among rings of size 5 or 6 in the molecule, we build a graph of rings that are “fused” 
     (sharing at least 2 atoms). 
  2. For each connected component we search for a subset of exactly 4 fused rings that:
       - contain exactly one 5-membered ring (the candidate D-ring) and three 6-membered rings;
       - have a union of atoms with between 15 and 19 carbon atoms (the classical steroid nucleus).
  3. In the candidate D-ring we then require that at least one carbon (with defined stereochemistry)
     is directly bonded in a single bond to an oxygen atom that appears to be a hydroxyl group 
     (i.e. the oxygen has at least one attached hydrogen).
     
NOTE: This is a heuristic and will not perfectly catch all steroids nor reject all non–steroids.
"""
from rdkit import Chem
from rdkit.Chem import rdchem
from itertools import combinations

def is_16beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 16β-hydroxy steroid based on its SMILES string using 
    a heuristic approach.

    A valid 16β–hydroxy steroid (heuristically) must:
      1. Possess a fused ring system extracted from rings of size 5 or 6 that contains exactly 4 rings,
         with 3 six–membered and 1 five–membered rings.
      2. This fused system (the steroid nucleus) should have roughly 15–19 carbon atoms.
      3. In the unique five–membered ring (candidate D–ring) there must be at least one carbon atom
         with defined stereochemistry (chiral tag not “unspecified”) that is directly bound via a single bond
         to an oxygen atom that appears to be an –OH group (the oxygen has at least one attached hydrogen).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule qualifies as a 16β-hydroxy steroid (by our heuristic), False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get all ring info; restrict to rings of size 5 or 6.
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()  # Each is a tuple of atom indices.
    rings_5_6 = []
    for ring in all_rings:
        if len(ring) in (5, 6):
            rings_5_6.append(set(ring))
    if not rings_5_6:
        return False, "No rings of size 5 or 6 found; cannot be a steroid nucleus"

    # Build a graph among rings (nodes are indices in rings_5_6). Two rings are fused if they share at least 2 atoms.
    num_rings = len(rings_5_6)
    ring_graph = {i: set() for i in range(num_rings)}
    for i in range(num_rings):
        for j in range(i+1, num_rings):
            if len(rings_5_6[i].intersection(rings_5_6[j])) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)

    # Find connected components among these rings.
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

    # Now, for each connected component with at least 4 rings, look for a combination of exactly 4 rings
    # that have the steroid nucleus features.
    candidate_nucleus = None
    candidate_D_ring = None
    for comp in components:
        if len(comp) < 4:
            continue
        comp_list = list(comp)
        # Try all combinations of 4 rings from this component.
        for comb in combinations(comp_list, 4):
            rings_subset = [rings_5_6[i] for i in comb]
            # Count ring sizes: exactly one 5-membered and three 6-membered.
            count5 = sum(1 for ring in rings_subset if len(ring)==5)
            count6 = sum(1 for ring in rings_subset if len(ring)==6)
            if count5 != 1 or count6 != 3:
                continue
            # Check connectivity within this subset:
            # Build a mini-graph for these 4 rings; each must be connected (share >=2 atoms) to at least one other.
            mini_connected = True
            for i in range(4):
                connected = False
                for j in range(4):
                    if i==j:
                        continue
                    if len(rings_subset[i].intersection(rings_subset[j])) >= 2:
                        connected = True
                        break
                if not connected:
                    mini_connected = False
                    break
            if not mini_connected:
                continue

            # Get union of atom indices in these 4 rings.
            union_atoms = set()
            for ring in rings_subset:
                union_atoms.update(ring)
            # Count carbons in the union.
            c_count = sum(1 for idx in union_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            # Classical steroid nucleus is around 17 carbons; allow a range from 15 to 19.
            if c_count < 15 or c_count > 19:
                continue

            # Identify the candidate D-ring (the only 5-membered ring in the subset)
            for ring in rings_subset:
                if len(ring) == 5:
                    # Also require that it is fused with at least one six–membered ring by (at least) 2 atoms.
                    has_fusion_with_6 = False
                    for ring2 in rings_subset:
                        if ring2 is ring or len(ring2) != 6:
                            continue
                        if len(ring.intersection(ring2)) >= 2:
                            has_fusion_with_6 = True
                            break
                    if not has_fusion_with_6:
                        continue
                    candidate_D_ring = ring
                    candidate_nucleus = union_atoms
                    break
            if candidate_nucleus is not None and candidate_D_ring is not None:
                break
        if candidate_nucleus is not None:
            break

    if candidate_nucleus is None:
        return False, "No fused 4–ring steroid nucleus (3 six–membered + 1 five–membered with ~17 carbons) found"

    # Now, within the candidate D–ring, check for a carbon with defined stereochemistry (chiral tag set)
    # that is directly bonded (via a single bond) to an oxygen that is likely –OH (has at least one hydrogen).
    found_16beta_oh = False
    reason_details = ""
    for idx in candidate_D_ring:
        atom = mol.GetAtomWithIdx(idx)
        # Select only carbon atoms
        if atom.GetAtomicNum() != 6:
            continue
        # Check if stereochemistry is defined.
        if atom.GetChiralTag() == rdchem.ChiralType.CHI_UNSPECIFIED:
            continue
        # Examine neighbors for an oxygen bound by a single bond.
        for nb in atom.GetNeighbors():
            if nb.GetAtomicNum() != 8:
                continue
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx())
            if bond is None or bond.GetBondType() != rdchem.BondType.SINGLE:
                continue
            # Check if the oxygen appears to be an –OH group;
            # typically this means it has at least one explicit or implicit hydrogen.
            if nb.GetTotalNumHs() < 1:
                continue
            # If we get here, we have a chiral carbon in the D–ring that bears an –OH.
            found_16beta_oh = True
            break
        if found_16beta_oh:
            break

    if not found_16beta_oh:
        return False, "Steroid nucleus found but no chiral carbon with an attached –OH in the 5–membered ring (D–ring) was detected"

    return True, ("Molecule contains a fused steroid nucleus (3 six–membered and 1 five–membered rings, ~17 carbons) "
                  "and the candidate D–ring has a chiral carbon bearing an –OH group (heuristic for 16β–hydroxy steroid)")

# Example usage:
if __name__ == '__main__':
    # Test with one true positive example: 16beta-hydroxytestosterone
    test_smiles = "C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1C[C@H](O)[C@@H]2O"
    result, reason = is_16beta_hydroxy_steroid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)