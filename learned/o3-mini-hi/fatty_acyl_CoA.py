"""
Classifies: CHEBI:37554 fatty acyl-CoA
"""
"""
Classifies: fatty acyl-CoA
Definition: An acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any fatty acid.
The molecule must contain a thioester group (i.e. R–C(=O)–S–) joined to a CoA fragment.
We use a minimal CoA fingerprint by looking for a substructure common to many acyl-CoA SMILES strings.
In addition, the acyl (fatty acid) part (the R-group attached to the carbonyl carbon) should be a moderately long (aliphatic) chain.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_fatty_acyl_CoA(smiles: str):
    """
    Determines if the given SMILES string corresponds to a fatty acyl-CoA.
    
    The criteria used are:
      1. Molecule must be valid.
      2. Must contain a thioester functional group (C(=O)S).
      3. The sulfur of the thioester must be linked to a coenzyme A (CoA) fragment.
         Here we approximate the CoA fragment by a minimal SMARTS ("SCCNC(=O)CCNC(=O)")
         found in many acyl-CoA molecules.
      4. The acyl part—the carbon chain attached to the carbonyl carbon (not the S side)—
         must be aliphatic and of a minimal length (we require at least 3 contiguous carbon atoms).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple where the Boolean indicates if SMILES is classified as a fatty acyl-CoA,
                     and the string gives a reason for the decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS for the thioester group: a carbonyl (C=O) directly bonded to a sulfur.
    thioester_smarts = "C(=O)S"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester (acyl-S) functional group found"
    
    # Define an approximate minimal SMARTS for the CoA moiety.
    # We look for the fragment: -SCCNC(=O)CCNC(=O)-
    # (Note: There is much variation in CoA SMILES but many acyl-CoA molecules contain this motif.)
    coa_smarts = "SCCNC(=O)CCNC(=O)"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    coa_hits = mol.GetSubstructMatches(coa_pattern)
    if not coa_hits:
        return False, "No CoA-like fragment found in molecule"
    # Flatten indices for easier look-up (all atoms that are part of a CoA hit)
    coa_atoms = set()
    for hit in coa_hits:
        for idx in hit:
            coa_atoms.add(idx)
    
    # Look for a thioester group match.
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group detected in molecule"
    
    # Define a helper function to perform a DFS over contiguous carbon atoms.
    def count_aliphatic_chain(start_idx, visited):
        count = 0
        stack = [start_idx]
        while stack:
            curr = stack.pop()
            if curr in visited:
                continue
            visited.add(curr)
            atom = mol.GetAtomWithIdx(curr)
            # We require that the atom is carbon and not part of the CoA fragment.
            if atom.GetAtomicNum() != 6 or curr in coa_atoms:
                continue
            count += 1
            for neighbor in atom.GetNeighbors():
                # Only follow bonds between carbons (either single or double)
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                    stack.append(neighbor.GetIdx())
        return count

    # Flag to record if an acceptable thioester link has been found.
    fatty_acyl_found = False
    reason_details = []
    
    # Loop over each thioester match and test if it bridges a fatty acyl and the CoA moiety.
    for match in thioester_matches:
        # In the pattern "C(=O)S", match[0] is the carbonyl carbon and match[1] is the sulfur.
        carbonyl_idx = match[0]
        sulfur_idx = match[1]
        # The sulfur should be part of the CoA fragment.
        if sulfur_idx not in coa_atoms:
            reason_details.append("Thioester sulfur not attached to a CoA fragment")
            continue
        # Get the carbonyl atom.
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # Identify the neighbor of the carbonyl that is NOT the sulfur and NOT the carbonyl oxygen.
        acyl_start = None
        for nbr in carbonyl_atom.GetNeighbors():
            # Skip the sulfur (our thioester connection) and highly electronegative atoms (O)
            if nbr.GetIdx() == sulfur_idx:
                continue
            if nbr.GetAtomicNum() == 8:  # likely the oxygen of the C=O group
                continue
            # We expect the fatty acyl chain to be carbon based.
            if nbr.GetAtomicNum() == 6:
                acyl_start = nbr.GetIdx()
                break
        if acyl_start is None:
            reason_details.append("No acyl chain found on the carbonyl side")
            continue
        
        # Count contiguous aliphatic (carbon) chain atoms starting at the acyl_start.
        visited = set()
        chain_length = count_aliphatic_chain(acyl_start, visited)
        if chain_length < 3:
            reason_details.append(f"Acyl chain too short (found chain length {chain_length})")
            continue
        # If we get here then we have found a thioester whose S is attached to CoA and the acyl chain is sufficiently long.
        fatty_acyl_found = True
        break

    if not fatty_acyl_found:
        # Return a reason that summarizes one or more failures
        if reason_details:
            return False, "; ".join(reason_details)
        else:
            return False, "No valid fatty acyl thioester fragment linked to a CoA fragment found"
    
    return True, "Molecule contains a thioester linking an acyl chain (fatty acid) with a CoA fragment"

# For testing purposes (you can uncomment the lines below to try some examples):
# example_smiles = "S(C(=O)CCC(CCCC(C)C)C)CCNC(=O)CCNC(=O)[C@H](O)C(COP(OP(OC[C@H]1O[C@@H](N2C3=NC=NC(N)=C3N=C2)C(O)[C@H]1OP(O)(O)=O)(O)=O)(O)=O)(C)C"
# result, reason = is_fatty_acyl_CoA(example_smiles)
# print(result, reason)