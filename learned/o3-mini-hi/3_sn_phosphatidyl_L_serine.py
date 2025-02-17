"""
Classifies: CHEBI:11750 3-sn-phosphatidyl-L-serine
"""
"""
Classifies: 3-sn-phosphatidyl-L-serine
Definition:
  A 3-sn-glycerophosphoserine compound having acyl substituents at the 1- and 2-hydroxy positions.
Improved heuristic:
    1. Parse the SMILES.
    2. Verify that the molecule contains at least one phosphorus atom.
    3. Check that a serine fragment [C](N)C(=O)O (ignoring stereochemistry) is present.
    4. Look for ester bonds (a carbonyl carbon that has a double-bonded O and a single-bonded O)
       and consider only those for which the –O– (linker) is “anchored” (reachable within 2 bonds) to 
       a phosphorus atom (the phosphoglycerol backbone).
    5. For each such ester bond (each candidate acyl chain), follow the acyl branch (starting from the 
       carbonyl carbon’s other carbon neighbor) restricting the traversal to carbons only; only count 
       branches with chain length ≥ 6 (chain credits the carbonyl carbon).
    6. Only if exactly 2 valid acyl chains (sn-1 and sn-2) are found is the molecule classified as a 
       3-sn-phosphatidyl-L-serine.
       
Note: This heuristic is not perfect. It may still misclassify molecules, but it reduces the false positives
by “anchoring” the acyl groups to the phosphate backbone.
"""

from rdkit import Chem

def is_3_sn_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a 3-sn-phosphatidyl-L-serine based on its SMILES string.
    The molecule must contain:
      (a) a phosphorus atom (for phosphate),
      (b) a serine fragment defined by [C](N)C(=O)O (ignoring chirality),
      (c) exactly two ester bonds that are attached (within two bonds) to a phosphorus atom,
          with the branch (acyl chain) coming off the carbonyl having a chain length of ≥6 carbons.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a 3-sn-phosphatidyl-L-serine, otherwise False.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that at least one phosphorus is present (for the phosphate group)
    if not any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Missing phosphorus (P) required for the phosphoglycerol head-group"
    
    # Check that a serine fragment is present.
    # We use a simplified SMARTS for serine: a carbon bound to an amino group and a carboxylic acid.
    serine_smarts = "C(N)C(=O)O"
    serine_pat = Chem.MolFromSmarts(serine_smarts)
    if serine_pat is None:
        return False, "Internal error processing serine SMARTS pattern"
    if not mol.HasSubstructMatch(serine_pat):
        return False, "Serine fragment (C(N)C(=O)O) not found in molecule"

    # Helper: Depth-first search to compute the longest contiguous carbon chain.
    # Starts from a given carbon atom (by index) and traverses only carbon atoms,
    # avoiding atoms in the banned set.
    def dfs_chain_length(start_idx, banned):
        max_length = 0
        stack = [(start_idx, 0, set())]  # (current_atom_idx, current_length, visited)
        while stack:
            cur_idx, cur_len, visited = stack.pop()
            if cur_idx in visited:
                continue
            visited = visited | {cur_idx}
            cur_len += 1
            if cur_len > max_length:
                max_length = cur_len
            atom = mol.GetAtomWithIdx(cur_idx)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr.GetAtomicNum() == 6 and nbr_idx not in banned:
                    stack.append((nbr_idx, cur_len, visited))
        return max_length

    # Helper: Check if an atom (typically an O atom in an ester bond) is connected (within 2 bonds) to a phosphorus.
    def is_connected_to_phosphorus(atom, max_depth=2):
        from collections import deque
        seen = set()
        queue = deque([(atom, 0)])
        while queue:
            current, depth = queue.popleft()
            if depth > max_depth:
                continue
            if current.GetAtomicNum() == 15:
                return True
            for nbr in current.GetNeighbors():
                if nbr.GetIdx() not in seen:
                    seen.add(nbr.GetIdx())
                    queue.append((nbr, depth+1))
        return False

    valid_acyl_count = 0
    counted_esters = set()  # to avoid double-counting the same ester motif
    
    # Now scan all bonds looking for ester motifs.
    # In an ester, a carbon (the carbonyl carbon) is double-bonded to one oxygen and single-bonded to another oxygen.
    for bond in mol.GetBonds():
        # Consider only single bonds.
        if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
            continue
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # We want one atom to be carbon (candidate "carbonyl carbon") and the other to be oxygen (the ester oxygen).
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 8:
            carbonyl = a1
            oester = a2
        elif a2.GetAtomicNum() == 6 and a1.GetAtomicNum() == 8:
            carbonyl = a2
            oester = a1
        else:
            continue
        
        # Check that the candidate carbonyl carbon has a double-bonded oxygen.
        dbl_bonded_found = False
        for nb in carbonyl.GetNeighbors():
            if nb.GetAtomicNum() == 8:
                b = mol.GetBondBetweenAtoms(carbonyl.GetIdx(), nb.GetIdx())
                if b is not None and b.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    dbl_bonded_found = True
                    break
        if not dbl_bonded_found:
            continue  # not an ester carbonyl

        # To consider this ester linkage, the oxygen (oester) should be attached
        # to the phosphoglycerol backbone. We require that oester is (within 2 bonds) connected to a phosphorus.
        if not is_connected_to_phosphorus(oester):
            continue

        # Get the acyl chain branch: from the carbonyl carbon, aside from the oester neighbor and the double-bonded O,
        # there should be a carbon that is the start of the acyl chain.
        acyl_start = None
        for nbr in carbonyl.GetNeighbors():
            # Skip the oester atom
            if nbr.GetIdx() == oester.GetIdx():
                continue
            # Also skip any oxygen bonded via double bond (the carbonyl oxygen)
            b = mol.GetBondBetweenAtoms(carbonyl.GetIdx(), nbr.GetIdx())
            if nbr.GetAtomicNum() == 8 and b is not None and b.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                continue
            if nbr.GetAtomicNum() == 6:
                acyl_start = nbr
                break
        if acyl_start is None:
            continue
        
        # Use the unique bond ID (tuple of sorted atom indices for the carbonyl–oester bond) to avoid double counting.
        bond_id = tuple(sorted((carbonyl.GetIdx(), oester.GetIdx())))
        if bond_id in counted_esters:
            continue
        counted_esters.add(bond_id)
         
        # Compute the acyl chain length. We include the carbonyl carbon (starting point)
        # and then traverse only carbon atoms. We do not cross back into the backbone (so ban the oester neighbor).
        chain_length = 1 + dfs_chain_length(acyl_start.GetIdx(), banned={carbonyl.GetIdx()})
        if chain_length >= 6:
            valid_acyl_count += 1

    if valid_acyl_count != 2:
        return False, f"Expected 2 acyl groups with sufficiently long chains anchored to a phosphate backbone; found {valid_acyl_count}"
    
    return True, "Molecule contains a phosphoserine head-group and 2 acyl groups (chain length ≥6) anchored to a phosphoglycerol backbone"

# Example usage (for testing; can be commented out):
if __name__ == "__main__":
    # Test with one of the provided examples: 1,2-distearoyl-sn-glycero-3-phosphoserine
    test_smiles = "CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)OC[C@H](N)C(O)=O)OC(=O)CCCCCCCCCCCCCCCCC"
    result, reason = is_3_sn_phosphatidyl_L_serine(test_smiles)
    print(result, reason)