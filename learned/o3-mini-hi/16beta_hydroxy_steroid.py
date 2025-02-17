"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
"""
Classifies: CHEBI:16β–hydroxy steroid
Definition: A 16–hydroxy steroid in which the hydroxy group at position 16 has a beta–configuration.
Heuristic approach (improved):
  1. Identify fused ring systems among rings of size 5–6. Two rings are fused if they share at least 2 atoms.
  2. Require at least 4 fused rings.
  3. In that fused system, require that one of the five–membered rings (candidate D–ring) is fused with at least one six–membered ring.
  4. Also, check that the total number of carbon atoms participating in the fused system is roughly in the range expected for a steroid nucleus (around 17 carbons).
  5. Finally, search within that five–membered ring for a carbon atom that is chiral (stereochemistry specified) and that is directly bonded (single bond) to an oxygen (–OH).
NOTE: This is only a heuristic method and will be imperfect.
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_16beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 16β-hydroxy steroid based on its SMILES string.

    A valid 16β–hydroxy steroid (heuristically) must:
      1. Contain a steroid nucleus identified as a fused ring system (rings of size 5 or 6)
         comprising at least 4 rings.
      2. Have within that fused system a five-membered ring (candidate D ring) that is fused
         with at least one six-membered ring.
      3. The fused nucleus should have roughly the expected number (~17) of carbon atoms.
      4. Possess in that five-membered ring a carbon atom with defined stereochemistry that bears
         an -OH group through a single bond.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule qualifies as 16β-hydroxy steroid (by our heuristic).
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring info (we use only rings of size 5 or 6)
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    rings_5_6 = [set(r) for r in all_rings if len(r) in (5, 6)]
    if not rings_5_6:
        return False, "No rings of size 5 or 6 found; cannot be a steroid nucleus."
    
    # Build a graph for rings: two rings are fused if they share at least 2 atoms.
    ring_graph = {i: [] for i in range(len(rings_5_6))}
    for i in range(len(rings_5_6)):
        for j in range(i+1, len(rings_5_6)):
            if len(rings_5_6[i].intersection(rings_5_6[j])) >= 2:
                ring_graph[i].append(j)
                ring_graph[j].append(i)
    
    # Find connected components (each component is a candidate fused system)
    seen = set()
    components = []
    for i in ring_graph:
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
    
    # Evaluate each component to see if it qualifies as a steroid nucleus.
    # Heuristics: at least 4 rings AND the nucleus (union of rings) should have ~17 carbons.
    steroid_component = None
    for comp in components:
        if len(comp) < 4:
            continue
        # Get unique atom indices in this component
        comp_atom_indices = set()
        for idx in comp:
            comp_atom_indices.update(rings_5_6[idx])
        # Count carbon atoms in the fused system
        c_count = sum(1 for idx in comp_atom_indices if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        # Expect between 15 and 19 carbons in the steroid nucleus (heuristic)
        if c_count < 15 or c_count > 19:
            continue

        # Check if at least one five-membered ring in the component is fused with a six-membered ring.
        comp_rings = [rings_5_6[i] for i in comp]
        candidate_D_ring = None
        for ring in comp_rings:
            if len(ring) == 5:
                # Look for a neighboring six-membered ring (sharing at least two atoms)
                found_neighbor = False
                for other in comp_rings:
                    if other is ring:
                        continue
                    if len(other) == 6 and len(ring.intersection(other)) >= 2:
                        found_neighbor = True
                        break
                if not found_neighbor:
                    continue
                # Now check within this five-membered ring for a chiral carbon with an -OH group.
                for idx in ring:
                    atom = mol.GetAtomWithIdx(idx)
                    if atom.GetAtomicNum() != 6:
                        continue
                    # Require stereochemistry to be set.
                    if atom.GetChiralTag() == rdchem.ChiralType.CHI_UNSPECIFIED:
                        continue
                    # Look among neighbors for an oxygen (atomic num 8) bound via a single bond.
                    oh_found = False
                    for nb in atom.GetNeighbors():
                        if nb.GetAtomicNum() == 8:
                            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx())
                            if bond and bond.GetBondType() == rdchem.BondType.SINGLE:
                                oh_found = True
                                break
                    if oh_found:
                        candidate_D_ring = ring
                        break
            if candidate_D_ring is not None:
                break
        
        if candidate_D_ring is not None:
            steroid_component = comp
            break
    
    if steroid_component is None:
        # Decide on appropriate message based on what failed.
        # It could be that no fused ring system was found with proper size, or no
        # five-membered ring with the expected chiral -OH feature was detected.
        return False, "No steroid nucleus meeting the fusion and carbon-count criteria with an appropriate five-membered ring bearing a chiral -OH was detected."
    
    return True, "Molecule contains a fused steroid nucleus (with ~17 carbons) and a five-membered ring having a chiral carbon with an -OH group (heuristic for 16β–hydroxy steroid)."

# Example usage:
# test_smiles = "C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1C[C@H](O)[C@@H]2O"  # Candidate: 16β-hydroxytestosterone
# print(is_16beta_hydroxy_steroid(test_smiles))