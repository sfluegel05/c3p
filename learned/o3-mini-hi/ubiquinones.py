"""
Classifies: CHEBI:16389 ubiquinones
"""
"""
Classifies: Ubiquinones – natural benzoquinones derived from 2,3-dimethoxy-5-methylbenzoquinone.
A valid ubiquinone must contain a redox‐active quinone core (cyclic diene-1,4-dione, as identified by
the SMARTS "O=C1C=CC(=O)C=C1"). Then, substituents on the core are inspected:
  • At least one acceptable oxygen substituent must be present (typically methoxy –OCH3 or hydroxy –OH).
  • In addition, there must be at least one alkyl substituent normally attached by carbon.
      - For a very short branch (fewer than 5 carbons), we allow even if a heteroatom is present.
      - For a long branch (5 or more carbons) the branch (intended to be a polyprenyl chain) must be acyclic,
        contain only carbon atoms and at least one non‐aromatic C=C bond.
  • No substituents other than these types should be attached directly to the quinone core.
Note: This algorithm is “best‐effort” and may fail on edge cases.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ubiquinones(smiles: str):
    """
    Determines if a molecule is classified as a ubiquinone.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): (True, explanation) if classified as a ubiquinone; or (False, reason) otherwise.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for quinone core: cyclic diene-1,4-dione pattern.
    core_smarts = "O=C1C=CC(=O)C=C1"
    core_query = Chem.MolFromSmarts(core_smarts)
    core_matches = mol.GetSubstructMatches(core_query)
    if not core_matches:
        return False, "Quinone core (cyclic diene-1,4-dione) not found"
    # For simplicity, use the first match as the core.
    core_atom_idxs = set(core_matches[0])
    
    # Counters for substituents.
    oxy_sub_count = 0   # acceptable oxygen substituents (methoxy/hydroxy)
    alkyl_count = 0     # acceptable alkyl substituents (methyl or polyprenyl chain)
    polyprenyl_detail = None  # if a long branch is found
    
    # We'll mark substituent atoms as visited so that we don’t re‐traverse them
    visited_substituents = set()
    
    def is_methyl(atom):
        # A carbon having at least 3 implicit hydrogens is considered methyl.
        return (atom.GetSymbol() == "C" and atom.GetTotalNumHs() >= 3)
    
    # Revised DFS that will traverse a branch (starting from a neighbor atom not in the core)
    # Returns:
    #   frag: set of atom indices in the branch (excluding the core)
    #   n_carb: number of carbon atoms in the branch
    #   has_double: True if any nonaromatic double bond (C=C) is encountered.
    #   has_ring: True if any atom in the branch is in a ring.
    #   has_hetero: True if any non-carbon (e.g. O, N, etc.) is encountered.
    def dfs_branch(start_atom, visited):
        stack = [start_atom]
        frag = set()
        n_carb = 0
        has_double = False
        has_ring = False
        has_hetero = False
        while stack:
            a = stack.pop()
            if a.GetIdx() in visited:
                continue
            visited.add(a.GetIdx())
            frag.add(a.GetIdx())
            if a.GetSymbol() == "C":
                n_carb += 1
            else:
                has_hetero = True  # encountered a heteroatom
            if a.IsInRing():
                has_ring = True
            for bond in a.GetBonds():
                nb = bond.GetOtherAtom(a)
                # Do not go back into the core.
                if nb.GetIdx() in core_atom_idxs:
                    continue
                # Continue along the branch if not visited.
                if nb.GetIdx() not in visited:
                    stack.append(nb)
                # If bond is a nonaromatic double, record it.
                if bond.GetBondTypeAsDouble() == 1 and not bond.GetIsAromatic():
                    has_double = True
        return frag, n_carb, has_double, has_ring, has_hetero

    # Now inspect each neighbor of the core (each direct substituent).
    for idx in core_atom_idxs:
        core_atom = mol.GetAtomWithIdx(idx)
        for nbr in core_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in core_atom_idxs or nbr_idx in visited_substituents:
                continue
            
            # Case A: If neighbor is oxygen, count as an acceptable oxygen substituent.
            if nbr.GetSymbol() == "O":
                # Get non-core neighbors from this oxygen.
                non_core_nbrs = [a for a in nbr.GetNeighbors() if a.GetIdx() not in core_atom_idxs]
                # No neighbor implies -OH (or deprotonated -O–).
                if len(non_core_nbrs) == 0:
                    oxy_sub_count += 1
                    visited_substituents.add(nbr_idx)
                    continue
                # If exactly one neighbor, check if it is a methoxy.
                if len(non_core_nbrs) == 1:
                    sub_atom = non_core_nbrs[0]
                    if sub_atom.GetSymbol() == "C" and is_methyl(sub_atom):
                        oxy_sub_count += 1
                        visited_substituents.add(nbr_idx)
                        visited_substituents.add(sub_atom.GetIdx())
                        continue
                    else:
                        # Otherwise, if the oxygen has H(s) attached, count as hydroxy.
                        if nbr.GetTotalNumHs() >= 1:
                            oxy_sub_count += 1
                            visited_substituents.add(nbr_idx)
                            continue
                        else:
                            # Otherwise, be generous and accept it.
                            oxy_sub_count += 1
                            visited_substituents.add(nbr_idx)
                            continue
                else:
                    # If oxygen is attached to >1 non-core neighbor, that is unexpected.
                    return False, f"Found oxygen substituent with {len(non_core_nbrs)} non-core neighbors which is not allowed."
            # Case B: If neighbor is carbon, then it is part of an alkyl substituent.
            elif nbr.GetSymbol() == "C":
                frag, n_carb, has_double, has_ring, has_hetero = dfs_branch(nbr, visited_substituents)
                # For very short branches (fewer than 5 carbons), we allow even if there are heteroatoms or rings.
                if n_carb < 5:
                    alkyl_count += 1
                    # If the branch is very short and qualifies as methyl, we note it.
                    # (We do not impose further restrictions here.)
                else:
                    # For longer branches, we require that:
                    #    – the branch is acyclic (no ring),
                    #    – contains only carbon atoms (no heteroatoms),
                    #    – and has at least one nonaromatic double bond.
                    if has_ring:
                        return False, f"Long substituent branch attached at atom {nbr_idx} contains a ring, which is not allowed."
                    if has_hetero:
                        return False, f"Long substituent branch attached at atom {nbr_idx} contains non-carbon atoms, which is not allowed."
                    if not has_double:
                        return False, f"Found a long substituent chain with {n_carb} carbons but no C=C double bond."
                    polyprenyl_detail = f"Found a polyprenyl chain with {n_carb} carbons containing C=C double bonds."
                    alkyl_count += 1
            else:
                # Any substituent attached to the core that is not O or C is not allowed.
                return False, f"Found substituent {nbr.GetSymbol()} attached to the core, which is not allowed."
    
    # Enforce minimal substituent requirements:
    if oxy_sub_count < 1:
        return False, f"Found {oxy_sub_count} acceptable oxygen substituent(s) on the quinone core; need at least 1 (typically matching a 2,3-dimethoxy pattern)."
    if alkyl_count < 1:
        return False, "No acceptable alkyl substituent found on the quinone core."
    
    # As an optional further check, rule out molecules that are too small.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150:
        return False, f"Molecular weight ({mol_wt:.1f} Da) too low for a typical ubiquinone."
    
    explanation = (f"Quinone core detected (matching {core_smarts}). Found {oxy_sub_count} acceptable oxygen substituent(s) "
                   f"and {alkyl_count} alkyl substituent(s). ")
    if polyprenyl_detail:
        explanation += polyprenyl_detail
    else:
        explanation += "No long polyprenyl side chain detected (this is acceptable for some ubiquinone homologues)."
    
    return True, explanation

# Example usage for testing:
if __name__ == "__main__":
    test_smiles = [
        "COC1=C(OC)C(=O)C(C\\C=C(/C)CCC=C(C)C)=C(C)C1=O",  # ubiquinone-2
        "C1(=C(C(C(=C(C1=O)OC)[O-])=O)C)C/C=C(/CC/C=C(/CC/C=C(/CC/C=C(/CC/C=C(/CCC=C(C)C)\\C)\\C)\\C)\\C)\\C",  # 3-demethyl-ubiquinone-7(1-)
        "C(C\\C(=C\\CC\\C(=C\\CC\\C(=C\\CC\\C(=C\\CC=1C(C(=C([O-])C(C1C)=O)OC)=O)\\C)\\C)\\C)/C=C(/CCC=C(C)C)\\C",  # 3-demethylubiquinone-6(1-)
        "COC1=C(OC)C(=O)C(C)=CC1=O",  # ubiquinone-0
        "COC1=C(OC)C(=O)C(CC=C(C)C)=C(C)C1=O",  # ubiquinone-1
        "COC1=C(OC)C(=O)C(C\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C)=C(C)C1=O",  # ubiquinone-8
        "O=C1C(OC)=C(OC)C(=O)C=C1[C@H](C(=O)C)C",  # 2,3-dimethoxy-5-(3-oxobutan-2-yl)cyclohexa-2,5-diene-1,4-dione (valid ubiquinone)
        "O(C=1C(=O)C([C@@H](C2=CC=C(O)C=C2)C=C)=CC(=O)C1OC)C",  # (R)-4'-Hydroxy-3,4-dimethoxydalbergione (valid ubiquinone)
        "O(C=1C(=O)C([C@@H](C2=CC=CC=C2)C=C)=CC(=O)C1OC)C",  # (R)-3,4-Dimethoxydalbergione (valid ubiquinone)
    ]
    
    for s in test_smiles:
        valid, reason = is_ubiquinones(s)
        print(f"SMILES: {s}\nResult: {valid}\nExplanation: {reason}\n")