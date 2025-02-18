"""
Classifies: CHEBI:16389 ubiquinones
"""
#!/usr/bin/env python
"""
Classifies: Ubiquinones – natural benzoquinones derived from 2,3-dimethoxy-5-methylbenzoquinone,
often with a polyprenyl side chain attached. Examples: ubiquinone-0, ubiquinone-2, 3-demethyl-ubiquinone-7, etc.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ubiquinones(smiles: str):
    """
    Determines if a molecule is a ubiquinone based on its SMILES string.
    Ubiquinones are defined as benzoquinones derived from 2,3-dimethoxy-5-methylbenzoquinone.
    They have a quinone ring (a cyclohexadiene-1,4-dione), with at least two methoxy substituents
    attached to the ring and usually a methyl group. Often, one substituent is a long polyprenyl chain.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a ubiquinone, False otherwise.
        str: Explanation/reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for a 1,4-benzoquinone core.
    # The SMARTS below is a simplified representation of a cyclohexadiene-1,4-dione.
    quinone_core = Chem.MolFromSmarts("O=C1C=CC(=O)C=C1")
    core_matches = mol.GetSubstructMatches(quinone_core)
    if not core_matches:
        return False, "Quinone core (cyclohexadiene-1,4-dione) not found"
    
    # For simplicity, take the first matching core
    core_indices = set(core_matches[0])
    
    # Initialize counts for substituents on the core.
    methoxy_count = 0  # count of -OCH3 substituents attached directly to the ring
    methyl_count = 0   # count of direct CH3 groups (not part of a methoxy)
    
    # To remember which neighbor atoms have been identified as substituents
    identified_substituents = set()
    
    # Loop over each atom in the quinone core to find substituents.
    for idx in core_indices:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in core_indices:
                continue  # skip atoms in the core
            # Check if the neighbor is an oxygen and could be part of a methoxy group.
            if nbr.GetSymbol() == "O":
                # We expect a methoxy to have the O attached to a CH3 (and only that non-core neighbor).
                nbr_neighbors = [a for a in nbr.GetNeighbors() if a.GetIdx() != idx]
                if len(nbr_neighbors) == 1:
                    attached = nbr_neighbors[0]
                    if attached.GetSymbol() == "C":
                        # use a simple check: if the attached carbon has three hydrogens (either implicit or explicit)
                        # then consider it a methyl group of a methoxy.
                        if attached.GetTotalNumHs() >= 3:
                            methoxy_count += 1
                            identified_substituents.add(nbr_idx)
                            identified_substituents.add(attached.GetIdx())
            # Otherwise check if the neighbor is a carbon that might be a direct methyl (CH3) substituent.
            elif nbr.GetSymbol() == "C":
                # Check if this carbon has three hydrogens and is only connected to the ring.
                # (Note: In longer chains the carbon will be connected to >1 heavy atom.)
                if nbr.GetTotalNumHs() >= 3 and nbr.GetDegree() == 1:
                    methyl_count += 1
                    identified_substituents.add(nbr_idx)
                    
    # We require at least two methoxy groups and at least one methyl substituent.
    if methoxy_count < 2:
        return False, f"Found only {methoxy_count} methoxy substituent(s); need at least 2"
    if methyl_count < 1:
        return False, "No methyl substituent directly attached to the core found"
    
    # Look for the presence of a side chain that might be a polyprenyl chain.
    # We search for substituents (neighbors of core atoms not already classified as methoxy or methyl)
    polyprenyl_found = False
    polyprenyl_reason = ""
    visited_substituents = set()
    
    def dfs_chain(start_atom, visited):
        """
        Depth-first search from a given atom (outside the core) to collect the connected substituent fragment.
        Returns a set of atom indices in that fragment and a boolean flag indicating whether a C=C (double bond)
        is present.
        """
        stack = [start_atom]
        frag_atoms = set()
        double_bond = False
        while stack:
            a = stack.pop()
            if a.GetIdx() in visited:
                continue
            visited.add(a.GetIdx())
            frag_atoms.add(a.GetIdx())
            for bond in a.GetBonds():
                nb = bond.GetOtherAtom(a)
                # Only follow if not part of the core.
                if nb.GetIdx() in core_indices:
                    continue
                if nb.GetIdx() not in visited:
                    stack.append(nb)
                if bond.GetBondTypeAsDouble() == 1:  # if bond type is DOUBLE
                    double_bond = True
        return frag_atoms, double_bond

    # Check each neighbor of a core atom that wasn't already classified as methoxy or methyl.
    for idx in core_indices:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in core_indices or nbr_idx in identified_substituents or nbr_idx in visited_substituents:
                continue
            # We require that the substituent is a carbon chain (at least one carbon)
            if nbr.GetSymbol() != "C":
                continue
            frag, has_dbl = dfs_chain(nbr, visited_substituents)
            # Count number of carbons in the fragment.
            n_carbons = sum(1 for i in frag if mol.GetAtomWithIdx(i).GetSymbol() == "C")
            # For a chain to be considered polyprenyl, require at least 5 carbon atoms and that a C=C double bond exists.
            if n_carbons >= 5 and has_dbl:
                polyprenyl_found = True
                polyprenyl_reason = f"Found substituent chain with {n_carbons} carbons and double bonds (indicative of a polyprenyl side chain)"
                break
        if polyprenyl_found:
            break

    # Build the explanation string.
    reason = f"Quinone core found with {methoxy_count} methoxy substituents and {methyl_count} methyl substituent(s). "
    if polyprenyl_found:
        reason += polyprenyl_reason
    else:
        reason += "No clear polyprenyl chain detected (side chain may be short or absent in some homologues)."
    
    return True, reason

# Example usage:
if __name__ == "__main__":
    # Test with a few example SMILES strings for ubiquinones.
    test_smiles = [
        "COC1=C(OC)C(=O)C(C\\C=C(/C)CCC=C(C)C)=C(C)C1=O",  # ubiquinone-2
        "COC1=C(OC)C(=O)C(C)=O",  # ubiquinone-0 (simplified)
    ]
    for s in test_smiles:
        valid, explanation = is_ubiquinones(s)
        print(f"SMILES: {s}\nResult: {valid}\nExplanation: {explanation}\n")