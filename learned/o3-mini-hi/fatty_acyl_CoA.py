"""
Classifies: CHEBI:37554 fatty acyl-CoA
"""
"""
Classifies: fatty acyl-CoA
Definition: An acyl-CoA results from the condensation of the thiol group of coenzyme A with the carboxy group of a fatty acid.
The molecule must contain a thioester group (i.e. C(=O)S) linking an acyl chain (fatty acid) on the carbonyl side to a CoA fragment.
We approximate a CoA fragment by requiring that the thioester sulfur is not isolated but is connected within a short bond distance
to a phosphorus atom (P) which is common in the phosphate groups of coenzyme A.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from collections import deque

def is_fatty_acyl_CoA(smiles: str):
    """
    Determines if the given SMILES string corresponds to a fatty acyl-CoA.
    
    The criteria used are:
      1. Molecule must be valid.
      2. Must contain a thioester group (C(=O)S) present in the molecule.
      3. The sulfur (S) from the thioester is expected to be linked (within a reasonable number of bonds)
         to at least one phosphorus (P) atom (common in the CoA moiety).
      4. An acyl chain must be present on the carbonyl side of the thioester. We require at least 3 contiguous 
         aliphatic carbon atoms starting from the carbon attached to the acyl chain.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple where the Boolean indicates if the molecule is classified as a fatty acyl-CoA,
                     and the string gives the reason for the decision.
    """
    # Parse the SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the thioester group (C(=O)S).
    thioester_smarts = "C(=O)S"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester (C(=O)S) functional group found"
    
    # Get indices for all phosphorus atoms in the molecule.
    phosphorus_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not phosphorus_indices:
        return False, "No phosphorus (P) atoms found; CoA moiety likely missing"
    
    # Helper function: breadth-first search to find if from a given atom index we can reach any phosphorus atom
    # within a maximum number of bonds (threshold).
    def is_connected_to_P(start_idx, threshold=15):
        visited = set([start_idx])
        queue = deque([(start_idx, 0)])
        while queue:
            curr, dist = queue.popleft()
            if dist > threshold:
                continue
            if curr in phosphorus_indices:
                return True, dist
            for nbr in mol.GetAtomWithIdx(curr).GetNeighbors():
                nidx = nbr.GetIdx()
                if nidx not in visited:
                    visited.add(nidx)
                    queue.append((nidx, dist + 1))
        return False, None

    # Helper function: count contiguous aliphatic carbon atoms starting from a given atom index.
    def count_aliphatic_chain(start_idx, excluded_idxs):
        count = 0
        visited = set()
        stack = [start_idx]
        while stack:
            curr = stack.pop()
            if curr in visited or curr in excluded_idxs:
                continue
            visited.add(curr)
            atom = mol.GetAtomWithIdx(curr)
            if atom.GetAtomicNum() != 6:  # must be carbon
                continue
            count += 1
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                    stack.append(nbr.GetIdx())
        return count

    reason_details = []
    fatty_acyl_found = False

    # Loop over each thioester match.
    for match in thioester_matches:
        # In our pattern "C(=O)S", match[0] is the carbonyl carbon (C) and match[1] is the sulfur (S).
        carbonyl_idx = match[0]
        sulfur_idx = match[1]
        
        # Check that the thioester sulfur is connected (within threshold bonds) to at least one phosphorus.
        connected, dist = is_connected_to_P(sulfur_idx, threshold=15)
        if not connected:
            reason_details.append("Thioester sulfur not connected to a phosphorus (CoA fragment) within 15 bonds")
            continue

        # Determine the acyl chain: from the carbonyl carbon, get the neighbor that is not the sulfur
        # and not the carbonyl oxygen.
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        acyl_start = None
        for nbr in carbonyl_atom.GetNeighbors():
            # Skip the sulfur atom of the thioester.
            if nbr.GetIdx() == sulfur_idx:
                continue
            # Skip oxygens, which are likely the carbonyl oxygen.
            if nbr.GetAtomicNum() == 8:
                continue
            if nbr.GetAtomicNum() == 6:
                acyl_start = nbr.GetIdx()
                break
        if acyl_start is None:
            reason_details.append("No acyl chain found on the carbonyl side")
            continue
        
        # Count contiguous carbons in the acyl chain (excluding atoms already involved in the thioester fragment).
        chain_length = count_aliphatic_chain(acyl_start, excluded_idxs={sulfur_idx, carbonyl_idx})
        if chain_length < 3:
            reason_details.append(f"Acyl chain too short (chain length {chain_length}); need at least 3 carbon atoms")
            continue

        # If we have a thioester with its sulfur properly linked to a phosphorus (CoA) and a sufficiently long acyl chain:
        fatty_acyl_found = True
        break

    if not fatty_acyl_found:
        if reason_details:
            return False, "; ".join(reason_details)
        else:
            return False, "No valid fatty acyl-CoA structure found"
    
    return True, "Molecule contains a thioester linking a fatty acyl chain with a CoA fragment (P-containing)"

# For testing purposes, you can uncomment the lines below:
# example_smiles = "S(C(=O)CCC(CCCC(C)C)C)CCNC(=O)CCNC(=O)[C@H](O)C(COP(OP(OC[C@H]1O[C@@H](N2C3=NC=NC(N)=C3N=C2)C(O)[C@H]1OP(O)(O)=O)(O)=O)(O)=O)(C)C"  # Dimethylnonanoyl-CoA
# result, reason = is_fatty_acyl_CoA(example_smiles)
# print(result, reason)