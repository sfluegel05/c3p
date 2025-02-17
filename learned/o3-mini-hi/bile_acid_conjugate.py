"""
Classifies: CHEBI:36249 bile acid conjugate
"""
"""
Classifies: Bile acid conjugates.
A bile acid conjugate is defined as a bile acid (a steroid nucleus made of 4 fused rings -
three six-membered rings and one five-membered ring) with a carboxylate side chain converted 
to an acyclic ester or amide and conjugated to a small functional fragment (e.g., glycine, taurine, 
sulfate, glucuronate, or sugar) that increases hydrophilicity.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdmolops

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid conjugate based on its SMILES string.
    
    Steps:
      1. Parse the SMILES string.
      2. Identify a fused steroid nucleus: a connected set of 4 rings (at least one 5-membered and three 
         six-membered rings) by grouping rings that share at least two atoms.
      3. Find acyclic ester (C(=O)O) or amide (C(=O)N) bonds whose carbonyl carbon is not in a ring 
         and is attached within the nucleus.
      4. Break the molecule at that bond to yield two fragments and then require that the small fragment 
         (the conjugate part) has a limited heavy atom count and matches a known conjugate SMARTS pattern.
    
    Returns:
        (bool, str): True and an explanation if classified as a bile acid conjugate; otherwise False with a reason.
    """
    
    def find_steroid_nucleus(mol):
        """
        Attempts to find a fused ring system that matches a bile acid nucleus.
        We check for a connected set of rings (rings are fused if they share at least 2 atoms)
        and then see if one connected group has exactly 4 rings with at least one 5-membered ring
        and three 6-membered rings.
        
        Returns:
            (set(atom indices), message) if found; otherwise (None, error message).
        """
        ring_info = mol.GetRingInfo()
        rings = ring_info.AtomRings()  # tuple of tuples of atom indices
        if len(rings) < 4:
            return None, f"Only {len(rings)} rings present; expect at least 4 for a bile acid nucleus."
        
        # Build graph among rings: nodes are ring indices; add edge if two rings share at least 2 atoms.
        n = len(rings)
        adj = {i: set() for i in range(n)}
        for i in range(n):
            for j in range(i+1, n):
                if len(set(rings[i]).intersection(rings[j])) >= 2:
                    adj[i].add(j)
                    adj[j].add(i)
        
        # Find connected components among rings (by indices).
        seen = set()
        components = []
        for i in range(n):
            if i in seen:
                continue
            comp = set()
            stack = [i]
            while stack:
                current = stack.pop()
                if current in comp:
                    continue
                comp.add(current)
                for neigh in adj[current]:
                    if neigh not in comp:
                        stack.append(neigh)
            seen |= comp
            components.append(comp)
        
        # Look for a component with exactly 4 rings.
        for comp in components:
            if len(comp) == 4:
                comp_rings = [rings[i] for i in comp]
                sizes = [len(ring) for ring in comp_rings]
                count5 = sum(1 for s in sizes if s == 5)
                count6 = sum(1 for s in sizes if s == 6)
                if count5 >= 1 and count6 >= 3:
                    atom_set = set()
                    for r in comp_rings:
                        atom_set |= set(r)
                    return atom_set, ("Found a proper bile acid nucleus: 4 fused rings "
                                      f"({count5} ring(s) with 5 atoms and {count6} ring(s) with 6 atoms).")
        return None, "No fused ring system matching bile acid nucleus criteria was found."
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Step 1: Identify the steroid nucleus.
    nucleus, nuc_msg = find_steroid_nucleus(mol)
    if nucleus is None:
        return False, f"Steroid nucleus not identified: {nuc_msg}"
    
    # Step 2: Define candidate linkage patterns: ester (C(=O)O) and amide (C(=O)N)
    linkage_patterns = [("ester", "C(=O)O"), ("amide", "C(=O)N")]
    
    # Pre-compile SMARTS patterns for known conjugate fragments:
    conjugate_patterns = [
        ("glycine", "[NX3;H2,H1]CC(=O)[O-]?"),      # glycine (amide-bound CH2COOH)
        ("taurine", "NCCS(=O)(=O)[O-]?"),            # taurine fragment
        ("sulfate", "[OS](=O)(=O)[O-]?"),             # sulfate group
        ("glucuronate", "OC1C(=O)[CH]([OH])[CH]([OH])[CH]([OH])[CH]1O"),  # simplified glucuronate
        ("sugar", "O[C@H]1OC(O)C(O)C(O)C(O)C1O")      # simple hexose
    ]
    conjugate_smarts = [(name, Chem.MolFromSmarts(smarts)) for name, smarts in conjugate_patterns]
    
    candidate_found = False
    candidate_reason = ""
    
    # Step 3: Scan for candidate linkage bonds.
    for link_name, link_pattern in linkage_patterns:
        patt = Chem.MolFromSmarts(link_pattern)
        if patt is None:
            continue
        matches = mol.GetSubstructMatches(patt)
        for match in matches:
            # Expect match has at least three atoms: the carbonyl carbon, its oxygen, and the conjugate heteroatom
            if len(match) < 3:
                continue
            carbonyl_idx = match[0]
            hetero_idx = match[2]
            
            # The carbonyl carbon should not be in a ring (i.e. it comes from the acyclic side-chain)
            if mol.GetAtomWithIdx(carbonyl_idx).IsInRing():
                continue
            # Further, require that the carbonyl carbon is part of the nucleus.
            if carbonyl_idx not in nucleus:
                continue
            # Also at least one neighbor (apart from the carbonyl oxygen) should be in the nucleus.
            nb_in_nuc = False
            for neighbor in mol.GetAtomWithIdx(carbonyl_idx).GetNeighbors():
                if neighbor.GetIdx() == match[1]:
                    continue
                if neighbor.GetIdx() in nucleus:
                    nb_in_nuc = True
                    break
            if not nb_in_nuc:
                continue
            
            # Get the bond between the carbonyl and the heteroatom.
            bond = mol.GetBondBetweenAtoms(carbonyl_idx, hetero_idx)
            if bond is None:
                continue
            bond_idx = bond.GetIdx()
            
            # Step 4: Fragment the molecule at this bond.
            try:
                frags = rdmolops.FragmentOnBonds(mol, [bond_idx], addDummies=True)
            except Exception:
                continue
            frag_mols = Chem.GetMolFrags(frags, asMols=True, sanitizeFrags=True)
            if len(frag_mols) != 2:
                continue
            
            # Decide which fragment is the steroid (nucleus) and which is the conjugate.
            steroid_frag = None
            conjugate_frag = None
            for frag in frag_mols:
                atom_ids = {atom.GetIdx() for atom in frag.GetAtoms()}
                # If the fragment contains many of the nucleus atoms, call it the steroid part.
                if len(atom_ids.intersection(nucleus)) >= len(nucleus) // 2:
                    steroid_frag = frag
                else:
                    conjugate_frag = frag
            if steroid_frag is None or conjugate_frag is None:
                continue
            
            # Require that the conjugate fragment is relatively small (heuristic: <=15 heavy atoms)
            if conjugate_frag.GetNumHeavyAtoms() > 15:
                continue
            
            # Check if the conjugate fragment matches one of the known patterns.
            for name, smarts_frag in conjugate_smarts:
                if smarts_frag is None:
                    continue
                if conjugate_frag.HasSubstructMatch(smarts_frag):
                    candidate_found = True
                    candidate_reason = (f"Found an acyclic {link_name} linkage splitting off a fragment "
                                        f"matching the conjugate pattern ({name}).")
                    break
            if candidate_found:
                break
        if candidate_found:
            break
    
    if not candidate_found:
        return False, ("No proper conjugate linkage was identified. Either no acyclic ester/amide bond connecting "
                       "a steroid-like nucleus (4 fused rings: 3 six-membered and 1 five-membered) to a small conjugate "
                       "fragment (glycine, taurine, sulfate, glucuronate or sugar) was detected.")
    
    return True, ("Molecule qualifies as a bile acid conjugate. " + nuc_msg + " " + candidate_reason)

# Example usage (for testing purposes only):
if __name__ == "__main__":
    # Example: one of the provided bile acid conjugate SMILES (Serine conjugated chenodeoxycholic acid)
    test_smiles = "O[C@H]1[C@]2([C@]3([C@@]([C@](CC3)([C@@H](CCC(=O)NC(CO)C(O)=O)C)[H])(CC[C@@]2([C@@]4([C@](C1)(C[C@H](O)CC4)[H])C)[H])C)[H])[H]"
    result, reason = is_bile_acid_conjugate(test_smiles)
    print("Result:", result)
    print("Reason:", reason)