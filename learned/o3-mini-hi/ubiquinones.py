"""
Classifies: CHEBI:16389 ubiquinones
"""
#!/usr/bin/env python
"""
Classifies: Ubiquinones – natural benzoquinones derived from 2,3-dimethoxy-5-methylbenzoquinone.
They have a redox‐active quinone core (a cyclic diene-1,4-dione) that normally carries one or two methoxy
substituents and at least one small alkyl (or polyprenyl) side chain.
This implementation first uses a SMARTS pattern to locate a 1,4‐benzoquinone core.
Then the substituents on the core are inspected:
  • At least one methoxy group (–OCH3) must be attached.
  • Also, at least one other small alkyl group (seen as a methyl –CH3, or a longer substituent which if ≥5 carbons
    must have a C=C double bond to be considered “polyprenyl”).
  • No other (unexpected) substituents are allowed.
If any long chain (≥5 carbons) is present it is inspected for a C=C bond.
If any substituent is “off‐pattern”, the molecule is not accepted.
Note: This algorithm is only one “best‐effort” approach.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ubiquinones(smiles: str):
    """
    Determines if a molecule is classified as a ubiquinone.
    The molecule must have a quinone core (a six-membered cyclic diene-1,4-dione) that is derived
    from 2,3-dimethoxy-5-methylbenzoquinone. In our implementation this means that:
      - A substructure matching "O=C1C=CC(=O)C=C1" is found in the molecule.
      - In the part of the molecule attached to this core, there is at least one methoxy substituent
        (–OCH3 on a ring atom).
      - There is at least one additional small alkyl substituent (either a methyl –CH3 or a longer chain)
        attached to the core. If a substituent chain has 5 or more carbons it must contain at least one
        C=C double bond (i.e. be a polyprenyl chain).
      - No substituent on the quinone core is “foreign” (for example, carrying halogens or extra heteroatoms).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): True if classified as a ubiquinone and an explanation, otherwise False and a reason.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS for a 1,4-benzoquinone core (cyclic diene-1,4-dione).
    core_smarts = "O=C1C=CC(=O)C=C1"
    core_query = Chem.MolFromSmarts(core_smarts)
    core_matches = mol.GetSubstructMatches(core_query)
    if not core_matches:
        return False, "Quinone core (cyclic diene-1,4-dione) not found"

    # For further analysis, take the first match.
    core_atom_idxs = set(core_matches[0])
    
    # We will now inspect substituents attached to the core.
    # Allowed substituents: methoxy group (OCH3), a methyl group (CH3), or a longer alkyl chain that 
    # – if ≥5 carbons – must contain a C=C double bond.
    methoxy_count = 0
    alkyl_count = 0  # counts methyl or polyprenyl (non-methoxy) groups
    polyprenyl_found = False
    polyprenyl_reason = ""
    
    # For keeping track of substituent atoms already claimed.
    visited_substituents = set()
    
    # Helper: check if an atom qualifies as a methyl group.
    def is_methyl(atom):
        return (atom.GetSymbol() == "C" and atom.GetTotalNumHs() >= 3)
    
    # Helper: do a DFS from a substituent neighbor and gather the fragment (do not cross back into the core)
    def dfs_fragment(start_atom, visited):
        stack = [start_atom]
        frag = set()
        has_double = False
        while stack:
            a = stack.pop()
            if a.GetIdx() in visited:
                continue
            visited.add(a.GetIdx())
            frag.add(a.GetIdx())
            for bond in a.GetBonds():
                nb = bond.GetOtherAtom(a)
                # do not cross into the core:
                if nb.GetIdx() in core_atom_idxs:
                    continue
                # add neighbor if not yet visited
                if nb.GetIdx() not in visited:
                    stack.append(nb)
                # Check bond order (nonaromatic double bonds)
                if bond.GetBondTypeAsDouble() == 1:
                    has_double = True
        return frag, has_double

    # Now loop over all atoms in the core and inspect each neighbor.
    # We also collect any neighbor that is not in the core.
    for idx in core_atom_idxs:
        core_atom = mol.GetAtomWithIdx(idx)
        for nbr in core_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in core_atom_idxs or nbr_idx in visited_substituents:
                continue
            # We now inspect the branch attached to the core.
            # If the neighbor is oxygen, check if this is a methoxy group.
            if nbr.GetSymbol() == "O":
                # Look at the lone neighbor (other than the core) on this oxygen.
                nbr_neighbors = [a for a in nbr.GetNeighbors() if a.GetIdx() != idx]
                if len(nbr_neighbors) == 1:
                    attached = nbr_neighbors[0]
                    if is_methyl(attached):
                        methoxy_count += 1
                        visited_substituents.add(nbr_idx)
                        visited_substituents.add(attached.GetIdx())
                        continue
                # If not a simple methoxy, mark as unsupported.
                return False, f"Unsupported oxygen substituent on the core (not a simple methoxy group) at atom {nbr_idx}"
            # If the neighbor is carbon, we further examine the fragment.
            elif nbr.GetSymbol() == "C":
                # Do DFS along this substituent branch.
                frag, has_double = dfs_fragment(nbr, visited_substituents)
                # Count number of carbon atoms in the fragment.
                n_carbons = sum(1 for atom_idx in frag if mol.GetAtomWithIdx(atom_idx).GetSymbol() == "C")
                # If it is a single carbon and qualifies as methyl:
                if n_carbons == 1 and is_methyl(nbr):
                    alkyl_count += 1
                # If it is a longer chain:
                elif n_carbons >= 5:
                    if not has_double:
                        return False, f"Found a long substituent chain with {n_carbons} carbons but with no C=C double bonds (likely a simple alkyl chain rather than a polyprenyl side chain)."
                    else:
                        polyprenyl_found = True
                        polyprenyl_reason = f"Found a polyprenyl chain with {n_carbons} carbons containing C=C double bonds."
                        alkyl_count += 1
                else:
                    # For chains of size 2-4 carbons, we expect them to be just a methyl or ethyl-type group.
                    alkyl_count += 1
            else:
                # Any substituent other than O or C is not allowed.
                return False, f"Found substituent {nbr.GetSymbol()} attached to the core, which is not allowed in ubiquinones."

    # Now enforce minimal requirements.
    if methoxy_count < 1:
        return False, f"Only found {methoxy_count} methoxy substituent(s) on the quinone core; need at least 1."
    if alkyl_count < 1:
        return False, "No alkyl (methyl or polyprenyl) substituent found on the quinone core."
    
    # Build a summary explanation.
    explanation = (f"Quinone core detected (as substructure {core_smarts}). Found {methoxy_count} methoxy group(s) and "
                   f"{alkyl_count} alkyl substituent(s). ")
    if polyprenyl_found:
        explanation += polyprenyl_reason
    else:
        explanation += "No long polyprenyl side chain detected (this is acceptable for some ubiquinone homologues)."
    
    return True, explanation

# Example usage for testing:
if __name__ == "__main__":
    test_smiles = [
        "COC1=C(OC)C(=O)C(C\\C=C(/C)CCC=C(C)C)=C(C)C1=O",  # ubiquinone-2
        "C1(=C(C(C(=C(C1=O)OC)[O-])=O)C)C/C=C(/CC/C=C(/CC/C=C(/CC/C=C(/CC/C=C(/CC/C=C(/CCC=C(C)C)\\C)\\C)\\C)\\C)\\C)\\C",  # 3-demethyl-ubiquinone-7(1-)
        "COC1=C(OC)C(=O)C(C)=CC1=O",  # ubiquinone-0
        "COC1=C(OC)C(=O)C(CC=C(C)C)=C(C)C1=O",  # ubiquinone-1
        "COC1=C(OC)C(=O)C(C\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C)=C(C)C1=O",  # ubiquinone-8
        # Some of the others in the provided outcomes could be added for further testing.
    ]
    
    for s in test_smiles:
        valid, reason = is_ubiquinones(s)
        print(f"SMILES: {s}\nResult: {valid}\nExplanation: {reason}\n")