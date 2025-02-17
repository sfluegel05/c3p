"""
Classifies: CHEBI:11750 3-sn-phosphatidyl-L-serine
"""
"""
Classifies: 3-sn-phosphatidyl-L-serine
Definition:
  A 3-sn-glycerophosphoserine compound having acyl substituents at the 1- and 2-hydroxy positions.
Heuristic approach:
1. Parse the SMILES and verify that the molecule contains a phosphorus atom (for the phospho group)
   and a chiral serine fragment defined by "[C@H](N)C(O)=O".
2. To find fatty acyl chains we iterate over all carbon atoms that are part of a carbonyl group:
   these carbons must have one double-bonded oxygen and one single-bonded oxygen (the ester oxygen).
3. For each such candidate, we identify the unique carbon substituent on the carbonyl carbon (the acyl chain,
   not the oxygen that links to the glycerol head-group) and measure the chain length.
4. The chain length is determined by a depth-first search that traverses only carbon atoms (atomic number 6)
   and does not cross back into the rest of the molecule (for example via the carbonyl carbon) once it has branched.
5. Only acyl chains with chain length ≥6 are accepted.
6. Exactly two acyl chains are expected.
"""

from rdkit import Chem

def is_3_sn_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a 3-sn-phosphatidyl-L-serine based on its SMILES string.
    The molecule must contain a glycerophosphoserine head-group (demonstrated by the presence of P and
    a serine fragment) and exactly two acyl chains (attached via ester bonds) that are sufficiently long
    (chain length of at least 6 carbon atoms, counting the carbonyl carbon).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a 3-sn-phosphatidyl-L-serine, else False.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check required atoms: phosphorus indicates the phosphate.
    if not any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Missing phosphorus (P) required for phosphoserine head-group"
    # Check for nitrogen (serine amino group)
    if not any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms()):
        return False, "Missing nitrogen (N) required for serine head-group"
    
    # Check for a serine fragment.
    # We use a chiral serine SMARTS: a carbon with amino and carboxyl groups.
    serine_smarts = "[C@H](N)C(O)=O"
    serine_pat = Chem.MolFromSmarts(serine_smarts)
    if serine_pat is None:
        return False, "Internal error processing serine SMARTS pattern"
    if not mol.HasSubstructMatch(serine_pat):
        return False, "Serine fragment ([C@H](N)C(O)=O) not found in molecule"
    
    # Define a helper function to compute the longest contiguous carbon chain
    # starting from a given atom (by index) while preventing to traverse any non-carbon.
    # We also pass an exclusion set so that once we leave the acyl branch we don't cross back.
    def dfs_chain_length(start_idx, banned):
        max_length = 0
        stack = [(start_idx, 0, set())]  # (current_atom_idx, current_length, visited)
        while stack:
            cur_idx, cur_len, visited = stack.pop()
            # Count current carbon if not seen already in current path.
            if cur_idx in visited:
                continue
            visited = visited | {cur_idx}
            cur_len += 1
            if cur_len > max_length:
                max_length = cur_len
            atom = mol.GetAtomWithIdx(cur_idx)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                # Only follow bonds to carbon and avoid any banned indices.
                if nbr.GetAtomicNum() == 6 and nbr_idx not in banned:
                    stack.append((nbr_idx, cur_len, visited))
        return max_length

    # Look for acyl chains attached via an ester group.
    # We iterate over carbons that are part of a carbonyl group (C=O).
    acyl_carbonyl_ids = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        # Identify double-bonded oxygen neighbor and collect single-bonded oxygens.
        dbl_oxy = None
        single_oxys = []
        for b in atom.GetBonds():
            nbr = b.GetOtherAtom(atom)
            if nbr.GetAtomicNum() == 8:
                if b.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    dbl_oxy = nbr
                elif b.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    single_oxys.append(nbr)
        # For an ester carbonyl, we expect one double-bonded O and at least one single-bonded O.
        # (In a typical ester, there is exactly one of each.)
        if dbl_oxy is None or len(single_oxys) < 1:
            continue
        # Identify the acyl chain branch.
        # The carbonyl carbon should have one carbon neighbor attached by a single bond
        # that is not the one coming from the ester linkage (which goes to O).
        acyl_neighbor = None
        for nbr in atom.GetNeighbors():
            # Exclude oxygens.
            if nbr.GetAtomicNum() == 6:
                # In an ester, there should be only one carbon neighbor which is the start of the acyl chain.
                acyl_neighbor = nbr
                break
        if acyl_neighbor is None:
            continue

        # Now, to compute the acyl chain length we consider the path starting from the acyl branch.
        # We do not want to "turn back" into the carbonyl carbon (which might lead to the ester oxygen branch);
        # so we set the carbonyl carbon index as banned.
        chain_length = 1 + dfs_chain_length(acyl_neighbor.GetIdx(), banned={atom.GetIdx()})
        # Use a threshold of 6 carbons (including the carbonyl carbon) to consider it a valid fatty acyl chain.
        if chain_length >= 6:
            acyl_carbonyl_ids.add(atom.GetIdx())
    
    # We expect exactly 2 acyl chains (one on sn-1 and one on sn-2).
    if len(acyl_carbonyl_ids) != 2:
        return False, f"Expected 2 acyl groups with sufficiently long chains; found {len(acyl_carbonyl_ids)}"
    
    return True, "Molecule contains a phosphoserine head-group and 2 acyl groups with long carbon chains"

# Example usage (can be commented out or used for testing):
if __name__ == "__main__":
    test_smiles = "CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)OC[C@H](N)C(O)=O)OC(=O)CCCCCCCCCCCCCCCCC"
    result, reason = is_3_sn_phosphatidyl_L_serine(test_smiles)
    print(result, reason)