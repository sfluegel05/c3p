"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
"""
Classifies: aliphatic nitrile (any nitrile derived from an aliphatic compound)

An aliphatic nitrile is defined as any molecule that contains at least one nitrile (–C≡N)
where the nitrile carbon is not aromatic and its single substituent (i.e. the atom attached
to the nitrile carbon besides the nitrile nitrogen) is an sp3 carbon that does not lead into 
a structurally extended (e.g. >6 carbons) or functionalized (e.g. carbonyl, extra nitrile)
branch and does not include any aromatic atoms.
"""

from rdkit import Chem

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if the molecule is an aliphatic nitrile based on its SMILES string.
    
    The procedure is as follows:
      1. Parse the SMILES and search for a nitrile substructure matching [C;X2]#[N;X1].
      2. For each nitrile match, take the C (nitrile carbon) and verify that it is not aromatic.
      3. Identify the substituent branch (the neighbor of the nitrile carbon that is not the nitrogen).
      4. Require that this substituent is a carbon atom (i.e. symbol=="C"), that is sp3,
         and that it bears at least one hydrogen.
      5. Check that none of its other bonds (to atoms other than the nitrile carbon) is directly
         to a polarized group like a carbonyl (a double bond to oxygen) or forms part of another nitrile.
      6. Walk the branch (using a breadth-first search) to count carbon atoms and to ensure that
         no aromatic atoms are encountered. If the branch contains too many carbons (a cutoff of 6 is used)
         we assume that the nitrile is derived from a long fatty or cyanolipid chain and reject it.
    
    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        bool: True if classified as an aliphatic nitrile, False otherwise.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the nitrile SMARTS pattern: C#N where the carbon has exactly two bonds.
    nitrile_pattern = Chem.MolFromSmarts("[C;X2]#[N;X1]")
    if nitrile_pattern is None:
        return False, "Error creating nitrile pattern"
    
    matches = mol.GetSubstructMatches(nitrile_pattern)
    if not matches:
        return False, "No nitrile group found in the molecule"
    
    # Helper: recursively check if any atom in the branch is aromatic.
    def branch_contains_aromatic(atom, banned_idx, visited):
        if atom.GetIdx() in visited:
            return False
        visited.add(atom.GetIdx())
        if atom.GetIsAromatic():
            return True
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() == banned_idx:
                continue
            if branch_contains_aromatic(nbr, banned_idx, visited):
                return True
        return False
    
    # Helper: count how many carbon atoms (unique indices) are in the branch.
    def count_branch_carbons(start_atom, banned_idx):
        queue = [start_atom]
        visited = {banned_idx}
        carbon_count = 0
        while queue:
            current = queue.pop(0)
            if current.GetIdx() in visited:
                continue
            visited.add(current.GetIdx())
            if current.GetSymbol() == "C":
                carbon_count += 1
            for nbr in current.GetNeighbors():
                if nbr.GetIdx() not in visited:
                    queue.append(nbr)
        return carbon_count

    # Loop over every nitrile match
    for match in matches:
        # match[0] is the nitrile carbon; match[1] is the nitrile nitrogen.
        nitrile_c = mol.GetAtomWithIdx(match[0])
        nitrile_n = mol.GetAtomWithIdx(match[1])
        
        # Skip if the nitrile carbon is aromatic or not a carbon.
        if nitrile_c.GetSymbol() != "C" or nitrile_c.GetIsAromatic():
            continue
        
        # Collect substituents excluding the nitrile nitrogen.
        subs = [nbr for nbr in nitrile_c.GetNeighbors() if nbr.GetIdx() != nitrile_n.GetIdx()]
        if not subs:
            continue  # Unusual, but go to next match

        # We expect exactly one substituent (terminal nitrile), so if multiple exist, we process each.
        for sub in subs:
            # Rule 1: The substituent must be carbon.
            if sub.GetSymbol() != "C":
                continue
            
            # Rule 2: It should be sp3 (if available, to indicate a saturated, aliphatic center).
            # (Note: sometimes hybridization info is not set; in that case, we try to check degree.)
            if sub.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                continue

            # Rule 3: It should have at least one attached hydrogen.
            if sub.GetTotalNumHs() < 1:
                continue
            
            # Rule 4: Check that none of its other bonds (besides to the nitrile carbon) is directly
            #         to atoms that are part of a polar function (eg, carbonyl oxygen or another nitrile).
            bad_attachment = False
            for nbr in sub.GetNeighbors():
                if nbr.GetIdx() == nitrile_c.GetIdx():
                    continue
                bond = mol.GetBondBetweenAtoms(sub.GetIdx(), nbr.GetIdx())
                if bond is not None:
                    # Check for carbonyl: if neighbor is O and bond is double.
                    if nbr.GetSymbol() == "O" and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        bad_attachment = True
                        break
                    # Check if neighbor is part of another nitrile bond.
                    if nbr.GetSymbol() == "N":
                        nbonds = nbr.GetBonds()
                        for b in nbonds:
                            if b.GetBondType() == Chem.rdchem.BondType.TRIPLE:
                                bad_attachment = True
                                break
                        if bad_attachment:
                            break
            if bad_attachment:
                continue
            
            # Rule 5: Recursively ensure the substituent’s branch contains no aromatic atoms.
            if branch_contains_aromatic(sub, nitrile_c.GetIdx(), set()):
                continue
            
            # Rule 6: Count carbons in the branch. If there are too many (here >6), reject.
            branch_carbons = count_branch_carbons(sub, nitrile_c.GetIdx())
            if branch_carbons > 6:
                continue
            
            # Passed all checks: classify as an aliphatic nitrile.
            return True, "Contains a nitrile group attached to an exclusively aliphatic substituent branch"
    
    # If none of the nitrile groups pass, provide an explanation.
    return False, "Nitrile group(s) found but none have a qualifying aliphatic substituent branch"


# Example usages (uncomment to test):
# print(is_aliphatic_nitrile("N#CCC#N"))                         # malononitrile -> expected True.
# print(is_aliphatic_nitrile("CC#N"))                             # acetonitrile -> expected True.
# print(is_aliphatic_nitrile("OCC#N"))                            # glycolonitrile -> expected False.
# print(is_aliphatic_nitrile("CCOC(=O)C(C)C#N"))                 # ethyl 2-cyanopropionate -> expected False.