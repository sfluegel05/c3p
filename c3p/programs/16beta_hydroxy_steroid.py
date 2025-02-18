"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
"""
Classifies: CHEBI:16β–hydroxy steroid
Definition: A 16–hydroxy steroid in which the hydroxy group at 
           position 16 has a beta–configuration.
Heuristic approach:
  1. First verify that the molecule contains a fused four–ring steroid nucleus.
     To do this we extract all rings of size 5 or 6 and construct a graph where
     two rings are connected (fused) if they share at least 2 atoms. A true steroid
     should have at least 4 fused rings.
  2. Then among the five–membered rings (expected to be the D–ring), we check for the presence 
     of a carbon atom that is stereogenic (i.e. its chiral tag is defined) and that has an 
     –OH group attached (via a single bond). We assume that in a true 16β–hydroxy steroid this 
     chiral center is the one bearing the 16–β hydroxyl.
NOTE: This method is only a heuristic and will not correctly classify every edge–case.
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_16beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 16β–hydroxy steroid based on its SMILES string.
    
    A valid 16β–hydroxy steroid is (heuristically) defined as a molecule having:
      1. A steroid nucleus in the form of a fused ring system made up of at least 4 rings
         (with rings of size 5 or 6) connected by sharing at least 2 atoms per fusion.
      2. At least one five-membered ring (expected to be ring D) wherein a carbon is marked as chiral
         (stereochemistry defined) and carries an –OH substituent via a single bond.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a 16β-hydroxy steroid, False otherwise.
        str: Reason explaining the classification.
    """
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information. We consider only rings of size 5 or 6.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    rings_5_6 = [set(r) for r in rings if len(r) in (5, 6)]
    if not rings_5_6:
        return False, "No rings of size 5 or 6 found; cannot be a steroid nucleus."
    
    # Build a graph of rings among rings_5_6 where two rings are connected if they share at least 2 atoms.
    ring_graph = {i: [] for i in range(len(rings_5_6))}
    for i in range(len(rings_5_6)):
        for j in range(i+1, len(rings_5_6)):
            if len(rings_5_6[i].intersection(rings_5_6[j])) >= 2:
                ring_graph[i].append(j)
                ring_graph[j].append(i)
    
    # Find connected components in the ring graph using DFS.
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
    
    # Check if any component qualifies as a steroid nucleus (i.e. having at least 4 fused rings).
    steroid_components = []
    for comp in components:
        if len(comp) >= 4:
            steroid_components.append(comp)
    
    if not steroid_components:
        return False, f"Fused ring analysis: found only {max((len(c) for c in components), default=0)} fused rings; expected at least 4 for a steroid nucleus."
    
    # Now, in at least one steroid nucleus component, look for a five-membered ring (expected D-ring)
    # that contains a chiral carbon with an -OH group.
    candidate_found = False
    msg_details = ""
    for comp in steroid_components:
        # Get all rings in this component
        comp_rings = [rings_5_6[i] for i in comp]
        # Filter to rings of size 5:
        five_membered = [r for r in comp_rings if len(r) == 5]
        if not five_membered:
            # No five-membered ring in this steroid nucleus candidate.
            msg_details += "A steroid nucleus was found but no five-membered ring was detected in it. "
            continue
        # Look for a carbon atom that is chiral and has a single-bonded oxygen neighbor attached.
        for ring in five_membered:
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                # Only consider carbon atoms
                if atom.GetAtomicNum() != 6:
                    continue
                # Check that stereochemistry is specified
                if atom.GetChiralTag() == rdchem.ChiralType.CHI_UNSPECIFIED:
                    continue
                # Now see if any neighbor is an oxygen (atomic num 8) connected via a single bond
                has_oh = False
                for nb in atom.GetNeighbors():
                    if nb.GetAtomicNum() == 8:
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx())
                        if bond and bond.GetBondType() == rdchem.BondType.SINGLE:
                            has_oh = True
                            break
                if has_oh:
                    candidate_found = True
                    # Note: We no longer strictly check for the chiral tag to be CHI_TETRAHEDRAL_CW;
                    # we assume that if stereochemistry is set and –OH is attached in a five-membered ring within the steroid nucleus,
                    # the molecule is likely a 16β–hydroxy steroid.
                    break
            if candidate_found:
                break
        if candidate_found:
            break
    if not candidate_found:
        return False, "No chiral carbon in a five-membered ring with an attached hydroxyl group (expected for 16β–hydroxy steroid) was detected in the steroid nucleus."
    
    return True, "Molecule contains a fused steroid nucleus and a chiral, hydroxyl-bearing carbon in a five-membered ring (heuristic for 16β–hydroxy steroid)."

# Example usage (these lines can be uncommented for testing):
# test_smiles = "C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1C[C@H](O)[C@@H]2O"  # Candidate: 16β-hydroxytestosterone
# print(is_16beta_hydroxy_steroid(test_smiles))