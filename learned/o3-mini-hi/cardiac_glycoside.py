"""
Classifies: CHEBI:83970 cardiac glycoside
"""
"""
Classifies: Cardiac glycoside
Definition: Steroid lactones containing sugar residues that act on the contractile force
            of the cardiac muscles.

This program checks for the presence of:
  (1) A steroid nucleus: approximated as a fused system of at least 4 rings 
      of size 5 or 6 that share at least 2 atoms between neighbors.
  (2) A lactone ring: we search among five‐membered rings for one that contains exactly one oxygen
      and a carbon with an exocyclic double bond to oxygen (i.e. a carbonyl group).
  (3) At least one sugar residue: a ring of size 5 (furanose) with at least one oxygen or size 6 (pyranose)
      with at least two oxygens that is not “buried” in the steroid core – that is, it either shares at most one atom 
      with the steroid nucleus or is adjacent (by a single bond) to it.
Note: Detecting fused ring systems and sugar residues is inherently approximate.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule is a cardiac glycoside based on its SMILES string.
    
    This function checks for:
      (1) A fused steroid nucleus, defined as a connected cluster of at least 4 rings (size 5 or 6)
          sharing at least 2 atoms between rings.
      (2) A lactone ring: a five-membered ring having one oxygen and containing a carbon that has an exocyclic double-bond
          to an oxygen (i.e. a carbonyl).
      (3) At least one sugar residue: detected from rings of size 5 or 6 with enough oxygens and that are attached 
          (or nearly attached) to the steroid nucleus. The lactone ring is excluded.
          
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): True if classified as a cardiac glycoside, along with a reason; otherwise False and a reason.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # -------------------------
    # (1) Detect the fused steroid nucleus.
    # Get all rings (as tuples of atom indices)
    ring_info = mol.GetRingInfo().AtomRings()
    candidate_rings = []
    for ring in ring_info:
        if len(ring) in (5, 6):
            candidate_rings.append(set(ring))
    if len(candidate_rings) < 4:
        return False, "Less than 4 rings of size 5 or 6 found; no steroid nucleus detected"
    
    # Build connectivity graph among candidate rings: two rings are connected if their intersection has at least 2 atoms.
    ring_adj = {i: set() for i in range(len(candidate_rings))}
    for i in range(len(candidate_rings)):
        for j in range(i+1, len(candidate_rings)):
            if len(candidate_rings[i].intersection(candidate_rings[j])) >= 2:
                ring_adj[i].add(j)
                ring_adj[j].add(i)
    
    # Collect connected components:
    visited = set()
    def dfs(node, current):
        current.add(node)
        visited.add(node)
        for nbr in ring_adj[node]:
            if nbr not in visited:
                dfs(nbr, current)
        return current

    steroid_atoms = set()
    steroid_found = False
    for i in range(len(candidate_rings)):
        if i not in visited:
            component = dfs(i, set())
            if len(component) >= 4:
                steroid_found = True
                # Union of all atoms in this fused set is our steroid nucleus
                for idx in component:
                    steroid_atoms.update(candidate_rings[idx])
                break
    if not steroid_found:
        return False, "No fused tetracyclic steroid nucleus detected"
    
    # -------------------------
    # (2) Detect the lactone ring.
    lactone_found = False
    lactone_ring_atoms = set()
    # Iterate through rings of size 5 as potential lactones.
    for ring in ring_info:
        if len(ring) != 5:
            continue
        # Count oxygen atoms in the ring.
        oxy_in_ring = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
        # Consider as lactone if exactly one oxygen is in the ring.
        if len(oxy_in_ring) != 1:
            continue
        # Now check if one of the carbon atoms in the ring has a double bond to an oxygen atom that is not a ring member.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # We expect a carbon that is part of a carbonyl
            if atom.GetAtomicNum() != 6:
                continue
            for bond in atom.GetBonds():
                # Check for a double bond to oxygen
                if bond.GetBondTypeAsDouble() >= 2.0:
                    nbr = bond.GetOtherAtom(atom)
                    if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in ring:
                        lactone_found = True
                        lactone_ring_atoms = set(ring)
                        break
            if lactone_found:
                break
        if lactone_found:
            break
    # Fallback SMARTS if nothing found by ring-analysis.
    if not lactone_found:
        lactone_smarts = "O=C1C=CCO1"
        lactone_pat = Chem.MolFromSmarts(lactone_smarts)
        if mol.HasSubstructMatch(lactone_pat):
            lactone_found = True
            # We take the first match as the lactone ring.
            lactone_ring_atoms = set(mol.GetSubstructMatches(lactone_pat)[0])
    if not lactone_found:
        return False, "No lactone ring (furanone/butenolide) found"
    
    # -------------------------
    # (3) Detect at least one sugar residue.
    sugar_found = False
    # We examine all rings (from ring_info) and select those of size 5 or 6.
    # Exclude the lactone ring.
    for ring in ring_info:
        ring_set = set(ring)
        if lactone_ring_atoms and ring_set == lactone_ring_atoms:
            continue  # do not confuse the lactone for a sugar.
        if len(ring) == 6:
            # For a pyranose-like ring, require at least 2 oxygen atoms.
            o_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if o_count < 2:
                continue
        elif len(ring) == 5:
            # For a furanose-like ring, require at least 1 oxygen.
            o_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if o_count < 1:
                continue
        else:
            continue
        # Now, decide whether this ring is “attached” (or adjacent) to the steroid nucleus.
        overlap = len(ring_set.intersection(steroid_atoms))
        attached = False
        if overlap <= 1:
            # If no atom in common, check if any atom in the ring is directly bonded to an atom in the steroid nucleus.
            if overlap == 0:
                for idx in ring_set:
                    atom = mol.GetAtomWithIdx(idx)
                    for nbr in atom.GetNeighbors():
                        if nbr.GetIdx() in steroid_atoms:
                            attached = True
                            break
                    if attached:
                        break
            else:
                attached = True
        if attached:
            sugar_found = True
            break
    # Fallback: try a loose standard sugar SMARTS if none found.
    if not sugar_found:
        sugar_smarts_loose = "O1CC(O)CC(O)C1"
        sugar_pat = Chem.MolFromSmarts(sugar_smarts_loose)
        if mol.HasSubstructMatch(sugar_pat):
            sugar_found = True
    if not sugar_found:
        return False, "No sugar residue (pyranose/furanose-like ring) found"
    
    # Optionally, add a molecular weight filter (cardiac glycosides are typically > 300 Da).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for typical cardiac glycosides"
    
    return True, "Contains a fused steroid nucleus, a lactone ring, and at least one sugar residue"

# (Optional) Example usage:
if __name__ == '__main__':
    # Example SMILES strings (from some true positives)
    test_smiles = [
        "[H]C(=O)O[C@H]1C[C@]2(O)[C@]3([H])CC[C@]4([H])C[C@H](CC[C@]4(C)[C@@]3([H])CC[C@]2(C)[C@@]1([H])C1=CC(=O)OC1)",  # gitaloxin
        "O[C@@]12[C@@](C)([C@H](CC1)C=3COC(C3)=O)[C@H](O)C[C@]4([C@]2(CC[C@]5([C@]4(C)CC[C@@H](C5)O[C@]6(C[C@H](O)[C@@H]([C@H](O6)C)O[C@]7(C[C@H](O)[C@@H]([C@H](O7)C)O[C@]8(C[C@H](OC(C)=O)[C@@H]([C@H](O8)C)O[C@]9(O[C@H](CO)[C@H]([C@@H]([C@H]9O)O)O)[H])[H])[H])[H])[H])[H])[H]",  # lanatoside C
    ]
    for smi in test_smiles:
        res, reason = is_cardiac_glycoside(smi)
        print("Result:", res, "| Reason:", reason)