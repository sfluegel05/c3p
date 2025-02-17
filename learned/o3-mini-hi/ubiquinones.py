"""
Classifies: CHEBI:16389 ubiquinones
"""
#!/usr/bin/env python
"""
Classifies: Ubiquinones – natural benzoquinones derived from 2,3-dimethoxy-5-methylbenzoquinone.
They possess a quinone ring (a cyclohexadiene-1,4-dione core) that carries at least one methoxy
substituent and a methyl group on the core. In many homologues one substituent is a polyprenyl chain,
which is defined as a carbon chain of at least 5 carbons that contains at least one C=C double bond.
We also want to exclude false positives with long alkyl chains that lack characteristic unsaturation.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ubiquinones(smiles: str):
    """
    Determines if a molecule is classified as a ubiquinone.
    For a molecule to be accepted it must contain:
      - A 1,4-benzoquinone core (cyclohexadiene-1,4-dione) – extra substituents allowed.
      - At least one methoxy substituent directly attached to the core. (Most derivatives have two,
        but 3-demethyl versions show only one.)
      - At least one methyl (CH3) substituent directly attached to the core.
      - If a substituent chain (non-methoxy, non-methyl) is present on the core with 5 or more carbons,
        it must contain at least one C=C (double bond) to qualify as a polyprenyl chain.
    Args:
        smiles (str): SMILES string of the molecule.
    Returns:
        bool, str: True and a reason if classified as a ubiquinone, False and a reason otherwise.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern to capture a cyclohexadiene-1,4-dione (1,4-benzoquinone core).
    # This pattern is kept flexible so that extra substituents on the ring are tolerated.
    quinone_core = Chem.MolFromSmarts("O=C1[CH]C([CH])C(=O)[CH]C1")  # using [CH] to loosely specify carbons
    core_matches = mol.GetSubstructMatches(quinone_core)
    if not core_matches:
        return False, "Quinone core (cyclohexadiene-1,4-dione) not found"
    
    # Use the first found core match.
    core_indices = set(core_matches[0])
    
    methoxy_count = 0  # –OCH3 groups
    methyl_count = 0   # directly attached CH3 groups
    identified_substituents = set()  # indices we have used for these small groups
    
    # Helper: decide if an atom is a CH3 (methyl) fragment.
    def is_methyl(atom):
        # We require a carbon that has at least three hydrogen attachments.
        return (atom.GetSymbol() == "C" and atom.GetTotalNumHs() >= 3)
    
    # Iterate over core indices and examine their neighbors.
    for idx in core_indices:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in core_indices:
                continue
            # Check for methoxy: an oxygen attached to the ring
            if nbr.GetSymbol() == "O":
                # Look for the other neighbor (should be CH3)
                # Sometimes an oxygen might be charged (e.g., [O-]); we still check its connection.
                nbr_neighbors = [a for a in nbr.GetNeighbors() if a.GetIdx() != idx]
                if len(nbr_neighbors) == 1:
                    attached = nbr_neighbors[0]
                    if is_methyl(attached):
                        methoxy_count += 1
                        identified_substituents.add(nbr_idx)
                        identified_substituents.add(attached.GetIdx())
            # Otherwise if neighbor is carbon, it might be a directly attached methyl or the start of a chain.
            elif nbr.GetSymbol() == "C":
                # If it is already identified as part of a methoxy group, skip it.
                if nbr_idx in identified_substituents:
                    continue
                # If it qualifies as a simple methyl (CH3) directly attached to the core, count it.
                # We do not restrict by degree here since sometimes stereochemistry or explicit atoms
                # may give a degree slightly higher than 1.
                if is_methyl(nbr):
                    methyl_count += 1
                    identified_substituents.add(nbr_idx)
    
    # Now look for any substituent chain attached to the core that wasn’t already classified as a small group.
    polyprenyl_found = False
    polyprenyl_reason = ""
    visited_substitution = set()
    
    # Depth-first search to discover the full fragment (chain) attached to the core.
    def dfs_chain(start_atom, visited):
        stack = [start_atom]
        frag_atoms = set()
        has_double = False
        while stack:
            a = stack.pop()
            if a.GetIdx() in visited:
                continue
            visited.add(a.GetIdx())
            frag_atoms.add(a.GetIdx())
            for bond in a.GetBonds():
                nb = bond.GetOtherAtom(a)
                # Do not cross into the core.
                if nb.GetIdx() in core_indices:
                    continue
                # Continue traversing.
                if nb.GetIdx() not in visited:
                    stack.append(nb)
                # Check if the bond is a double bond.
                if bond.GetBondTypeAsDouble() == 1:
                    has_double = True
        return frag_atoms, has_double

    # Iterate over core neighbors (not already identified as small substituents) to see if there is a chain.
    chain_error = None
    for idx in core_indices:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in core_indices or nbr_idx in identified_substituents or nbr_idx in visited_substitution:
                continue
            # Consider only carbon substituents as potential chain starters.
            if nbr.GetSymbol() != "C":
                continue
            frag, has_dbl = dfs_chain(nbr, visited_substitution)
            # Count carbons in the fragment.
            n_carbons = sum(1 for i in frag if mol.GetAtomWithIdx(i).GetSymbol() == "C")
            # If the fragment is small, assume it is just a small substituent already counted.
            if n_carbons >= 5:
                # For a polyprenyl side chain, we require at least one double bond.
                if has_dbl:
                    polyprenyl_found = True
                    polyprenyl_reason = f"Found substituent chain with {n_carbons} carbons and C=C double bonds (indicative of a polyprenyl side chain)"
                else:
                    # If a long chain is present but without any double bonds, that is not allowed.
                    chain_error = f"Found long substituent chain with {n_carbons} carbons but no C=C double bonds (likely a non-isoprenoid alkyl chain)"
                # In either case, we break after the first chain found.
                break
        if chain_error or polyprenyl_found:
            break

    # In our classification, we require at least one methoxy group.
    # (Most ubiquinones have two but the de-methylated forms have one.)
    if methoxy_count < 1:
        return False, f"Found only {methoxy_count} methoxy substituent(s); need at least 1"
    if methyl_count < 1:
        return False, "No methyl substituent directly attached to the core found"
    
    # If a long chain is attached and it is not a polyprenyl chain, then reject the molecule.
    if chain_error:
        return False, chain_error

    # Build the explanation.
    reason = f"Quinone core detected with {methoxy_count} methoxy substituent(s) and {methyl_count} methyl substituent(s). "
    if polyprenyl_found:
        reason += polyprenyl_reason
    else:
        reason += "No polyprenyl side chain detected (side chain may be short or absent in some homologues)."
    
    return True, reason

# Example usage for testing:
if __name__ == "__main__":
    # List of test SMILES strings (names commented) from the provided outcomes.
    test_smiles = [
        # True positives:
        "COC1=C(OC)C(=O)C(C\\C=C(/C)CCC=C(C)C)=C(C)C1=O",  # ubiquinone-2
        "COC1=C(OC)C(=O)C(C)=CC1=O",  # ubiquinone-0
        "COC1=C(OC)C(=O)C(CC=C(C)C)=C(C)C1=O",  # ubiquinone-1
        "COC1=C(OC)C(=O)C(C\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C)=C(C)C1=O",  # ubiquinone-8
        "O(C=1C(=O)C(C\\C=C(\\CC\\C=C(\\CC\\C=C(\\CCC=C(C)C)/C)/C)/C)=C(C(=O)C1OC)C)C",  # Coenzyme Q4
        "COC1=C(OC)C(=O)C(C\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C)=C(C)C1=O",  # ubiquinone-5
        "COC1=C(OC)C(=O)C(C\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C)=C(C)C1=O",  # coenzyme Q10
        "O=C1C(OC)=C(OC)C(=O)C(=C1CC=C(CCCC(=O)C)C)",  # Pseudoalteromone A
        "O=C1C(OC)=C(OC)C(=O)C(=C1C/C=C(/CC/C=C(/C[C@H](O)C=C(C)C)/C)/C)",  # Antroquinonol N
        "COC1=C(OC)C(=O)C(C\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C)=C(C)C1=O",  # ubiquinone-7
        "COC1=C(OC)C(=O)C(C\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C)=C(C)C1=O",  # ubiquinone-9
        # False negatives (should be classified as ubiquinones):
        "C1(=C(C(C(=C(C1=O)OC)[O-])=O)C)C/C=C(/CC/C=C(/CC/C=C(/CC/C=C(/CC/C=C(/CC/C=C(/CCC=C(C)C)\\C)\\C)\\C)\\C)\\C)\\C",  # 3-demethyl-ubiquinone-7(1-)
        "C(C\\C(=C\\CC\\C(=C\\CC\\C(=C\\CC\\C(=C\\CC=1C(C(=C([O-])C(C1C)=O)OC)=O)\\C)\\C)\\C)\\C)/C=C(/CCC=C(C)C)\\C",  # 3-demethylubiquinone-6(1-)
        "C1(=C(C(C(=C(C1=O)OC)[O-])=O)C)C/C=C(/CC/C=C(/CC/C=C(/CC/C=C(/CC/C=C(/CC/C=C(/CCC=C(C)C)\\C)\\C)\\C)\\C)\\C)\\C",  # 3-demethylubiquinone-8(1-)
        "O=C1C(OC)=C(OC)C(=O)C=C1[C@H](C(=O)C)C",  # 2,3-dimethoxy-5-(3-oxobutan-2-yl)...
        "O(C=1C(=O)C([C@@H](C2=CC=C(O)C=C2)C=C)=CC(=O)C1OC)C",  # (R)-4'-Hydroxy-3,4-dimethoxydalbergione
        "O(C=1C(=O)C([C@@H](C2=CC=CC=C2)C=C)=CC(=O)C1OC)C",  # (R)-3,4-Dimethoxydalbergione
        "S(C=1C(=O)C(OC)=C(OC)C(C1CC(=O)OC)=O)C",  # Coptirhoquinone A
        "O(C=1C(=O)C(=C(CC=C(C)C)C(=O)C1OC)CC=C(C)C)C=2C(=C(O)C=C(O)C2)C(OC)=O",  # Atrovirinone
        "COC1=C(O)C(=O)C(C)=C(C/C=C(/C)CC/C=C(/C)CC/C=C(/C)CC/C=C(/C)CC/C=C(/C)CC/C=C(/C)CC/C=C(/C)CCC=C(C)C)C1=O",  # 3-demethylubiquinone-9
        "O=C1C(OC)=C(OC)C(=O)C=C1C(C(O)C)C",  # 5-(3-hydroxybutan-2-yl)-2,3-dimethoxy...
    ]
    
    for s in test_smiles:
        valid, explanation = is_ubiquinones(s)
        print(f"SMILES: {s}\nResult: {valid}\nExplanation: {explanation}\n")