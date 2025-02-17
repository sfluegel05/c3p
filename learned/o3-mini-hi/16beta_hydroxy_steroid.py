"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
"""
Classifies: CHEBI:16β–hydroxy steroid
Definition: A 16–hydroxy steroid in which the hydroxy group at position 16 has a beta–configuration.
Improved heuristic:
  1. Extract rings of size 5 or 6 and group fused rings via shared atoms.
  2. In each connected component, search for exactly 4 rings (3 six–membered and 1 five–membered)
     that are fused (each shares ≥2 atoms with at least one other) and whose union contains exactly 17 carbons.
  3. In the unique 5–membered (“D–ring”) from that nucleus, look for at least one carbon that:
       a. Has defined chirality,
       b. Has connectivity (when counting neighbors in the nucleus) consistent with a peripheral (non–junction) carbon,
       c. Is directly bound via a single bond to an oxygen atom that carries at least one hydrogen (–OH).
Note: This heuristic does not cover every steroid variant perfectly but is tuned toward reducing false positives.
"""
from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem import AllChem
from itertools import combinations

def is_16beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 16β–hydroxy steroid via an improved heuristic.
    The molecule must possess a fused four–ring system (steroid nucleus) composed of three six–membered and one five–membered rings.
    The union of atoms in these rings must yield exactly 17 carbon atoms.
    In the unique five–membered (D–ring) at least one carbon with defined stereochemistry (and non–junction connectivity)
    must be directly bonded via a single bond to an oxygen atom having at least one hydrogen (i.e. an –OH group).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule qualifies as a 16β–hydroxy steroid (by our heuristic), False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure stereochemistry is assigned.
    Chem.AssignStereochemistry(mol, force=True)
    
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()  # tuples of atom indices
    rings_5_6 = []   # consider only rings with 5 or 6 members
    for ring in all_rings:
        if len(ring) in (5, 6):
            rings_5_6.append(set(ring))
    if not rings_5_6:
        return False, "No rings of size 5 or 6 found; cannot be a steroid nucleus"
    
    # Build a graph over rings: nodes are rings and edges are drawn if the rings share at least 2 atoms.
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
    candidate_D_ring = None  # will be the only five-membered ring in the nucleus
    # Explore each connected component for a valid combination of 4 rings
    for comp in components:
        if len(comp) < 4:
            continue
        comp_list = list(comp)
        # Examine all 4-ring combinations in the component.
        for comb in combinations(comp_list, 4):
            rings_subset = [rings_5_6[i] for i in comb]
            # Must have exactly 1 five-membered ring and 3 six-membered rings.
            count5 = sum(1 for ring in rings_subset if len(ring) == 5)
            count6 = sum(1 for ring in rings_subset if len(ring) == 6)
            if count5 != 1 or count6 != 3:
                continue
            # They should be mutually fused; each ring must share >=2 atoms with at least one other.
            mini_connected = True
            for i in range(4):
                if not any(i != j and len(rings_subset[i].intersection(rings_subset[j])) >= 2 
                           for j in range(4)):
                    mini_connected = False
                    break
            if not mini_connected:
                continue
            # Get the union of all atoms in these rings.
            union_atoms = set()
            for ring in rings_subset:
                union_atoms.update(ring)
            # Count carbon atoms among these.
            c_count = sum(1 for idx in union_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            if c_count != 17:
                continue
            
            # Identify the candidate D–ring: the unique 5–membered ring.
            for ring in rings_subset:
                if len(ring) == 5:
                    # To be sure, this 5–membered ring should be fused with at least one six–membered ring.
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
        return False, "No fused 4–ring steroid nucleus (3 six–membered and 1 five–membered rings with exactly 17 carbons) found"
    
    # Now, within the candidate D–ring, look for a proper candidate chiral carbon bearing an –OH.
    found_16beta_oh = False
    for idx in candidate_D_ring:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() != 6:
            continue
        # Require that the chiral tag is defined.
        if atom.GetChiralTag() == rdchem.ChiralType.CHI_UNSPECIFIED:
            continue

        # Count the number of neighbors that are part of the nucleus.
        neighbors_in_nucleus = sum(1 for nb in atom.GetNeighbors() if nb.GetIdx() in candidate_nucleus)
        # In a classical steroid, the C16 carbon (in the D–ring) is not a full fusion point.
        # We require it has exactly one neighbor in the nucleus (the ring junction) since its other bond (C17) lies outside.
        if neighbors_in_nucleus != 1:
            continue

        # Search among neighbors for an oxygen that appears to be in an –OH.
        for nb in atom.GetNeighbors():
            if nb.GetAtomicNum() != 8:
                continue
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx())
            if bond is None or bond.GetBondType() != rdchem.BondType.SINGLE:
                continue
            # Require that the oxygen has at least one hydrogen explicitly attached.
            if nb.GetTotalNumHs(includeNeighbors=True) < 1:
                continue
            found_16beta_oh = True
            break
        if found_16beta_oh:
            break

    if not found_16beta_oh:
        return False, "Steroid nucleus found but no chiral carbon (with appropriate connectivity) in the D–ring bearing an –OH was detected"
    
    return True, ("Molecule contains a fused steroid nucleus (3 six–membered and 1 five–membered rings with exactly 17 carbons) "
                  "and a candidate D–ring carbon (with defined stereochemistry and proper connectivity) bearing an –OH group, "
                  "consistent with a 16β–hydroxy steroid.")

# Example usage:
if __name__ == '__main__':
    # Test example: 16beta-hydroxytestosterone
    test_smiles = "C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1C[C@H](O)[C@@H]2O"
    result, reason = is_16beta_hydroxy_steroid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)