"""
Classifies: CHEBI:83970 cardiac glycoside
"""
"""
Classifies: Cardiac glycoside
Definition: Steroid lactones containing sugar residues that act on the contractile force
            of the cardiac muscles.

Heuristic detection checks for:
  (1) A fused steroid nucleus. We examine rings of size 5 or 6 and then build a connectivity
      graph among rings sharing at least 2 atoms. A connected component of at least 4 such rings
      is taken as the steroid nucleus.
  (2) A lactone ring. We look among 5-membered rings for one that has exactly one oxygen and that 
      contains at least one carbon bearing an exocyclic (i.e. not in-ring) double-bond to an oxygen (a carbonyl).
      In addition we require that the lactone ring is fused to the steroid core.
  (3) At least one sugar residue. We search for rings of size 5 (furanose) or 6 (pyranose) that have 
      at least one oxygen atom. To avoid mis‐identifying the lactone as a sugar, any ring that exactly matches 
      the lactone ring is skipped. Also, the sugar ring should be “peripheral” – either sharing at most one atom 
      with the steroid nucleus or having a direct bond from one of its atoms to an atom in the steroid nucleus.
      
Note: These criteria are heuristic. There will inevitably be mis‐classifications.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule is a cardiac glycoside based on its SMILES string.
    
    Returns:
        (bool, str): A tuple containing True with a reason if the molecule is classified as a 
                     cardiac glycoside; otherwise False with an explanation.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # -------------------------
    # (1) Detect the fused steroid nucleus.
    # Get all ring systems (set of atom indices) from RDKit ring info.
    ring_info = mol.GetRingInfo().AtomRings()
    candidate_rings = []
    for ring in ring_info:
        if len(ring) in (5, 6):
            candidate_rings.append(set(ring))
    if len(candidate_rings) < 4:
        return False, "Less than 4 rings (size 5 or 6) found; no steroid nucleus detected"
    
    # Build connectivity among candidate rings: two rings are connected if they share at least 2 atoms.
    ring_adj = {i: set() for i in range(len(candidate_rings))}
    for i in range(len(candidate_rings)):
        for j in range(i + 1, len(candidate_rings)):
            if len(candidate_rings[i].intersection(candidate_rings[j])) >= 2:
                ring_adj[i].add(j)
                ring_adj[j].add(i)
                
    # Collect connected components via DFS; require one component has at least 4 rings.
    visited = set()
    def dfs(node, comp):
        comp.add(node)
        visited.add(node)
        for nbr in ring_adj[node]:
            if nbr not in visited:
                dfs(nbr, comp)
        return comp
    steroid_atoms = set()
    steroid_found = False
    for i in range(len(candidate_rings)):
        if i not in visited:
            comp = dfs(i, set())
            if len(comp) >= 4:
                steroid_found = True
                for idx in comp:
                    steroid_atoms.update(candidate_rings[idx])
                break
    if not steroid_found:
        return False, "No fused tetracyclic steroid nucleus detected"
    
    # -------------------------
    # (2) Detect the lactone ring.
    lactone_found = False
    lactone_ring_atoms = set()
    # Look among 5-membered rings.
    for ring in ring_info:
        if len(ring) != 5:
            continue
        # Count oxygen atoms in the ring.
        oxy_in_ring = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
        if len(oxy_in_ring) != 1:
            continue
        # Look for a carbon in the ring that has a double bond (exocyclic carbonyl) to an oxygen outside the ring.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue
            for bond in atom.GetBonds():
                # Check that the bond order is double and the neighbor is oxygen not in the ring.
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
    # Fallback: try a loose SMARTS match if nothing found by ring‐analysis.
    if not lactone_found:
        lactone_smarts = "O=C1C=CCO1"
        lactone_pat = Chem.MolFromSmarts(lactone_smarts)
        if mol.HasSubstructMatch(lactone_pat):
            lactone_found = True
            lactone_ring_atoms = set(mol.GetSubstructMatches(lactone_pat)[0])
    if not lactone_found:
        return False, "No lactone ring (furanone/butenolide) found"
    # Improve lactone quality: require that the lactone is fused to the steroid nucleus.
    connected_to_steroid = False
    for idx in lactone_ring_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in steroid_atoms:
                connected_to_steroid = True
                break
        if connected_to_steroid:
            break
    if not connected_to_steroid:
        return False, "Lactone ring found but not fused to steroid nucleus"
    
    # -------------------------
    # (3) Detect at least one sugar residue.
    sugar_found = False
    # Examine all rings (from ring_info) of size 5 or 6, excluding the lactone ring.
    for ring in ring_info:
        ring_set = set(ring)
        if lactone_ring_atoms and ring_set == lactone_ring_atoms:
            continue  # do not confuse the lactone with a sugar.
        if len(ring) not in (5, 6):
            continue
        # For many deoxysugars, require at least one oxygen (instead of 2).
        o_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        if o_count < 1:
            continue
        
        # Check if the sugar ring is attached "peripherally" to the steroid nucleus.
        # If it shares at most one atom with the steroid or if (when no atom is shared) any atom in the ring
        # is directly bonded to a steroid atom.
        overlap = len(ring_set.intersection(steroid_atoms))
        attached = False
        if overlap <= 1:
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
    # Fallback: search with a loose SMARTS pattern for a sugar (here pyranose-like).
    if not sugar_found:
        sugar_smarts_loose = "O1CC(O)CC(O)C1"
        sugar_pat = Chem.MolFromSmarts(sugar_smarts_loose)
        if mol.HasSubstructMatch(sugar_pat):
            sugar_found = True
    if not sugar_found:
        return False, "No sugar residue (pyranose/furanose-like ring) found"
    
    # Optionally, add a molecular weight filter; many cardiac glycosides are above ~300 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for typical cardiac glycosides"
    
    return True, "Contains a fused steroid nucleus, a lactone ring (fused to the nucleus), and at least one sugar residue"

# (Optional) Example usage:
if __name__ == '__main__':
    # A few example SMILES strings (these are taken from the provided positive and negative examples)
    test_smiles = [
        # True positives:
        "[H]C(=O)O[C@H]1C[C@]2(O)[C@]3([H])CC[C@]4([H])C[C@H](CC[C@]4(C)[C@@]3([H])CC[C@]2(C)[C@@]1([H])C1=CC(=O)OC1",  # gitaloxin
        "O[C@@]12[C@@](C)([C@H](CC1)C=3COC(C3)=O)[C@H](O)C[C@]4([C@]2(CC[C@]5([C@]4(C)CC[C@@H](C5)O[C@]6(C[C@H](O)[C@@H]([C@H](O6)C)O[C@]7(C[C@H](O)[C@@H]([C@H](O7)C)O[C@]8(C[C@H](OC(C)=O)[C@@H]([C@H](O8)C)O[C@]9(O[C@H](CO)[C@H]([C@@H]([C@H]9O)O)O)[H])[H])[H])[H])[H])[H])[H]",  # lanatoside C
        # False positive (should be rejected):
        "O=C1C2=C(C(=O)C[C@@H]3[C@@]2(CC[C@@H](C3(C)C)O)C)[C@@]4(C(=O)C[C@]5([C@]4(C1)C)O[C@]6(OC(=O)[C@@H](C6)C)CC5C)C",  # (17S,23S)-17,23-epoxy-...
    ]
    for smi in test_smiles:
        res, reason = is_cardiac_glycoside(smi)
        print("SMILES:", smi)
        print("Result:", res, "| Reason:", reason)
        print("----------------------------------------")