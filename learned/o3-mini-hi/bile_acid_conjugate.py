"""
Classifies: CHEBI:36249 bile acid conjugate
"""
#!/usr/bin/env python
"""
Classifies: Bile acid conjugates.
A bile acid conjugate is defined as a bile acid (a molecule with a steroid‐like, fused four–ring nucleus including at least one five–membered and one six–membered ring, and a carboxylate side chain that has been converted to an amide or ester)
that is conjugated via that linkage to a functional group that increases hydrophilicity (e.g. glycine, taurine, sulfate, glucuronate or a simple sugar).
This implementation uses heuristic SMARTS‐based filters plus an explicit check that the linkage bond “splits” the molecule into a steroid core and a small conjugate fragment.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdmolops

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid conjugate based on its SMILES string.
    
    The function performs the following steps:
      1. Parse the SMILES and check that the molecule has at least 4 rings,
         including at least one five–membered and at least one six–membered ring.
      2. Search for an acyclic ester or amide linkage (C(=O)O or C(=O)N) where
         the carbonyl carbon is not in a ring and is attached (outside the carbonyl oxygen) 
         to a ring – this is used as a proxy for the carboxylate side chain of a bile acid.
      3. For each candidate linkage, “break” the bond and isolate the fragment that is not
         part of the steroid nucleus. Then check whether the fragment matches one of the
         expected conjugate fragment patterns.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if classified as a bile acid conjugate, False otherwise.
        str: Reason explaining the classification decision.
    """
    # ---------- Step 1: Parse SMILES and check for steroid nucleus ---------- #
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Use ring info to get the rings present.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if len(atom_rings) < 4:
        return False, f"Insufficient rings: found {len(atom_rings)}, expected at least 4 for a bile acid nucleus."
    
    ring_sizes = [len(r) for r in atom_rings]
    if 5 not in ring_sizes:
        return False, "No five-membered ring found (expected one D-ring in a bile acid nucleus)."
    if 6 not in ring_sizes:
        return False, "No six-membered ring found (expected at least one six-membered ring in the steroid nucleus)."
    
    # ---------- Step 2: Find candidate conjugation linkages by searching for acyclic ester/amide bonds ---------- #
    # We will use SMARTS to find C(=O)O and C(=O)N linkages.
    ester_smarts = "C(=O)O"  # for esters
    amide_smarts = "C(=O)N"  # for amides
    patterns = [("ester", ester_smarts), ("amide", amide_smarts)]
    
    candidate_found = False
    candidate_reason = ""
    
    # Pre-compile conjugate fragment SMARTS patterns.
    # These patterns are simplified and ignore chiral specifications.
    conjugate_patterns = [
        ("glycine", "[NX3;H2,H1]CC(=O)[O-]?"),  # glycine-like fragment (amide-bound CH2COOH)
        ("taurine", "NCCS(=O)(=O)[O-]?"),         # taurine fragment
        ("sulfate", "[OS](=O)(=O)[O-]?"),          # sulfate group
        ("glucuronate", "OC1C(=O)[CH]([OH])[CH]([OH])[CH]([OH])[CH]1[OH]"),  # simplified glucuronate
        ("sugar", "OC[C@@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1O")  # a simple hexose pattern
    ]
    conjugate_smarts = [(name, Chem.MolFromSmarts(smarts)) 
                        for name, smarts in conjugate_patterns]
    
    # Loop over the two linkage patterns.
    for link_name, link_smarts in patterns:
        patt = Chem.MolFromSmarts(link_smarts)
        if patt is None:
            continue
        matches = mol.GetSubstructMatches(patt)
        for match in matches:
            # For both ester and amide, match[0] is the carbonyl carbon.
            carbonyl_idx = match[0]
            # Check that the carbonyl atom is NOT in a ring.
            if mol.GetAtomWithIdx(carbonyl_idx).IsInRing():
                continue
            
            # Identify the "conjugate side" atom.
            # In our SMARTS, match[2] is the heteroatom (O for ester, N for amide).
            conjugate_atom_idx = match[2]
            
            # Check that the carbonyl carbon is attached to some ring atom.
            carbon_neighbors = mol.GetAtomWithIdx(carbonyl_idx).GetNeighbors()
            attached_to_ring = False
            for nb in carbon_neighbors:
                # Skip the carbonyl oxygen (which is double-bonded).
                if nb.GetIdx() == match[1]:
                    continue
                if nb.IsInRing():
                    attached_to_ring = True
                    break
            if not attached_to_ring:
                continue  # likely not coming from a steroid side chain
            
            # Now that we have a candidate linkage, break just the bond between the carbonyl carbon and the conjugate side.
            bond = mol.GetBondBetweenAtoms(carbonyl_idx, conjugate_atom_idx)
            if bond is None:
                continue
            bond_idx = bond.GetIdx()
            
            # Fragment the molecule on that bond.
            try:
                frags = rdmolops.FragmentOnBonds(mol, [bond_idx], addDummies=True)
            except Exception as e:
                continue
            frag_mols = Chem.GetMolFrags(frags, asMols=True, sanitizeFrags=True)
            # We expect two fragments.
            if len(frag_mols) != 2:
                continue
            
            # Decide which fragment likely contains the steroid nucleus.
            # We assume the steroid nucleus will be the fragment having >= 4 rings.
            steroid_frag = None
            conjugate_frag = None
            for frag in frag_mols:
                ri = frag.GetRingInfo().AtomRings()
                if len(ri) >= 4:
                    steroid_frag = frag
                else:
                    conjugate_frag = frag
            if steroid_frag is None or conjugate_frag is None:
                continue  # not a clear split
            
            # Now check if the conjugate fragment contains one of the known conjugate patterns.
            for c_name, c_patt in conjugate_smarts:
                if c_patt is None:
                    continue
                if conjugate_frag.HasSubstructMatch(c_patt):
                    candidate_found = True
                    candidate_reason = (f"Found a candidate {link_name} linkage that splits off a fragment "
                                        f"matching a conjugate pattern ({c_name}).")
                    break
            if candidate_found:
                break
        if candidate_found:
            break

    if not candidate_found:
        return False, ("No proper conjugate linkage found. Either no acyclic ester/amide bond connecting "
                       "a steroid-like core to a small conjugate fragment (glycine, taurine, sulfate, glucuronate or sugar) was detected.")
    
    return True, ("Contains a steroid nucleus (>=4 rings with at least one five- and one six-membered ring) "
                  "and a conjugation linkage that yields a fragment matching a known conjugate group. " + candidate_reason)

# Example usage (for testing):
if __name__ == "__main__":
    # Test with one bile acid conjugate example SMILES
    # (Serine conjugated chenodeoxycholic acid):
    test_smiles = "O[C@H]1[C@]2([C@]3([C@@]([C@](CC3)([C@@H](CCC(=O)NC(CO)C(O)=O)C)[H])(CC[C@@]2([C@@]4([C@](C1)(C[C@H](O)CC4)[H])C)[H])C)[H])[H]"
    result, reason = is_bile_acid_conjugate(test_smiles)
    print("Result:", result)
    print("Reason:", reason)