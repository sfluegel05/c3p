"""
Classifies: CHEBI:36249 bile acid conjugate
"""
"""
Classifies: Bile acid conjugates.
A bile acid conjugate is defined here as a molecule that contains a steroid nucleus (4 fused rings:
3 six-membered rings and 1 five-membered ring) with a carboxylate side chain that has been
converted via an acyclic ester or amide linkage to a small conjugate fragment (such as glycine,
taurine, sulfate, glucuronate, a sugar, or other small polar functional groups).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdmolops

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid conjugate based on its SMILES string.
    
    The strategy is:
      1. Parse the SMILES.
      2. Identify the steroid nucleus as a fused ring system of 4 rings (3 six-membered and 1 five-membered).
      3. Search for an acyclic ester (C(=O)O) or amide (C(=O)N) bond. For a valid conjugation bond,
         we require that the carbonyl carbon is not in a ring and that it is connected (via at least one neighbor)
         to the steroid nucleus.
      4. Break the molecule at that bond and select the small fragment containing the conjugate.
         Verify that the small fragment (by heavy atom count) contains a known conjugate motif.
    
    Returns:
        (bool, str): True with explanation if the molecule qualifies, otherwise False and a reason.
    """
    
    def find_steroid_nucleus(mol):
        """
        Look for a fused ring system that is candidate for the bile acid nucleus.
        We require that one connected set of rings (shared by at least two atoms each) has exactly 4 rings,
        containing at least one 5-membered and three 6-membered rings.
        
        Returns:
            (set(atom indices), message) if found; otherwise (None, error message).
        """
        ring_info = mol.GetRingInfo()
        rings = ring_info.AtomRings()
        if len(rings) < 4:
            return None, f"Only {len(rings)} rings found; expecting 4 fused rings for a steroid nucleus."
        
        # Build connectivity between rings (they are 'fused' if they share at least 2 atoms).
        n = len(rings)
        adj = {i: set() for i in range(n)}
        for i in range(n):
            for j in range(i+1, n):
                if len(set(rings[i]).intersection(rings[j])) >= 2:
                    adj[i].add(j)
                    adj[j].add(i)
        
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
        
        # Look for a component that has exactly 4 rings and the proper ring-size distribution.
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
                    return atom_set, (f"Found bile acid nucleus: 4 fused rings ({count5} five-membered and {count6} six-membered).")
        return None, "No fused ring system matching bile acid nucleus criteria was found."
    
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Step 1: Identify the steroid nucleus.
    nucleus, nuc_msg = find_steroid_nucleus(mol)
    if nucleus is None:
        return False, f"Steroid nucleus not identified: {nuc_msg}"
    
    # Step 2: Define candidate linkage patterns: ester (C(=O)O) and amide (C(=O)N)
    linkage_patterns = [("ester", "C(=O)[O,N]")]  # using a more generic pattern
    # Pre-compile SMARTS for known conjugate fragments.
    # These patterns are simplified and may match a fragment having a polar small group.
    conjugate_patterns = [
        ("glycine", "[NX3;H2,H1]CC(=O)[O-]?"),
        ("taurine", "NCCS(=O)(=O)[O-]?"),
        ("sulfate", "[OS](=O)(=O)[O-]?"),
        ("glucuronate", "OC1C(=O)[CH]([OH])[CH]([OH])[CH]([OH])[CH]1O"),
        ("sugar", "O[C@H]1OC(O)C(O)C(O)C(O)C1O")
    ]
    conjugate_smarts = [(name, Chem.MolFromSmarts(smarts)) for name, smarts in conjugate_patterns]
    
    candidate_found = False
    candidate_reason = ""
    
    # Step 3: Find each acyclic ester or amide bond candidate.
    for link_name, link_smarts in linkage_patterns:
        link_mol = Chem.MolFromSmarts(link_smarts)
        if link_mol is None:
            continue
        matches = mol.GetSubstructMatches(link_mol)
        for match in matches:
            # Assume the first atom in the match is the carbonyl carbon,
            # the second atom is the carbonyl oxygen (or N) attached via double bond,
            # and the third is the heteroatom of the conjugated part.
            if len(match) < 3:
                continue
            carbonyl_idx = match[0]
            # We take the heteroatom on the other side of the bond. Sometimes the order may vary,
            # so examine the neighbors of the carbonyl that are not the carbonyl oxygen.
            carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
            if carbonyl_atom.IsInRing():
                continue  # skip if carbonyl carbon is in a ring
            
            hetero_idx = None
            for neigh in carbonyl_atom.GetNeighbors():
                # skip the carbonyl oxygen (or the double-bonded atom)
                if neigh.GetAtomicNum() in (8, 7):  # oxygen or nitrogen are candidates
                    # choose the neighbor which is not in the same fragment as the nucleus, but connected to it indirectly.
                    hetero_idx = neigh.GetIdx()
                    break
            if hetero_idx is None:
                continue
            
            # Check that the carbonyl atom has at least one neighbor (besides the heteroatom) that is within the steroid nucleus.
            connected_to_nucleus = False
            for neigh in carbonyl_atom.GetNeighbors():
                if neigh.GetIdx() == hetero_idx:
                    continue
                if neigh.GetIdx() in nucleus:
                    connected_to_nucleus = True
                    break
            if not connected_to_nucleus:
                continue  # not linking to the nucleus
            
            # Now we fragment on the bond between the carbonyl carbon and the heteroatom.
            bond = mol.GetBondBetweenAtoms(carbonyl_idx, hetero_idx)
            if bond is None or bond.IsInRing():
                continue  # want an acyclic bond
            bond_idx = bond.GetIdx()
            try:
                fragmented = rdmolops.FragmentOnBonds(mol, [bond_idx], addDummies=True)
            except Exception:
                continue
            frag_mols = Chem.GetMolFrags(fragmented, asMols=True, sanitizeFrags=True)
            if len(frag_mols) != 2:
                continue
            
            # Decide which fragment is likely the steroid nucleus: the fragment that shares more atoms with our nucleus.
            steroid_frag = None
            conjugate_frag = None
            for frag in frag_mols:
                frag_atom_ids = {atom.GetAtomMapNum() for atom in frag.GetAtoms() if atom.HasProp("_MolFileAtomStereo")}
                # If the fragment contains at least one atom from the nucleus, consider heavy atom count.
                frag_indices = {atom.GetIdx() for atom in frag.GetAtoms()}
                if len(frag_indices.intersection(nucleus)) > 0:
                    # Heuristic: larger fragment is likely the nucleus.
                    if steroid_frag is None or frag.GetNumHeavyAtoms() > steroid_frag.GetNumHeavyAtoms():
                        steroid_frag = frag
                    else:
                        conjugate_frag = frag
                else:
                    # The fragment that does not share with the nucleus is the conjugate.
                    conjugate_frag = frag
            
            if steroid_frag is None or conjugate_frag is None:
                continue
            
            # Ensure the conjugate fragment is relatively small.
            if conjugate_frag.GetNumHeavyAtoms() > 20:
                continue
            
            # Check the conjugate fragment against our known patterns.
            for name, smarts_frag in conjugate_smarts:
                if smarts_frag is None:
                    continue
                if conjugate_frag.HasSubstructMatch(smarts_frag):
                    candidate_found = True
                    candidate_reason = (f"Found an acyclic {link_name} linkage with a small fragment matching the conjugate pattern ({name}).")
                    break
            if candidate_found:
                break
        if candidate_found:
            break
    
    if not candidate_found:
        return False, ("No proper conjugate linkage was identified. The criteria require an acyclic ester/amide bond "
                       "linking a steroid nucleus (4 fused rings) to a small conjugate fragment (e.g., glycine, taurine, "
                       "sulfate, glucuronate or sugar) that increases hydrophilicity.")
    
    return True, ("Molecule qualifies as a bile acid conjugate. " + nuc_msg + " " + candidate_reason)

# Example usage (for testing purposes only):
if __name__ == "__main__":
    # Test with one of the provided bile acid conjugate SMILES:
    test_smiles = "O[C@H]1[C@]2([C@]3([C@@]([C@](CC3)([C@@H](CCC(=O)NC(CO)C(O)=O)C)[H])(CC[C@@]2([C@@]4([C@](C1)(C[C@H](O)CC4)[H])C)[H])C)[H])[H]"
    result, reason = is_bile_acid_conjugate(test_smiles)
    print("Result:", result)
    print("Reason:", reason)