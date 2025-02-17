"""
Classifies: CHEBI:37554 fatty acyl-CoA
"""
"""
Classifies: fatty acyl-CoA
Definition: An acyl-CoA results from the condensation of the thiol group of coenzyme A with the carboxy group of a fatty acid.
The molecule must contain a thioester group (C(=O)S) linking an acyl chain (fatty acid) on the carbonyl side to a CoA fragment.
We approximate the CoA fragment by requiring the presence both of at least one phosphorus (P) atom along with a purine (adenine) substructure.
"""

from rdkit import Chem
from collections import deque

def is_fatty_acyl_CoA(smiles: str):
    """
    Determines if the given SMILES string corresponds to a fatty acyl-CoA.
    Our criteria are:
      1. The SMILES is valid.
      2. Must contain a thioester group (C(=O)S).
      3. The thioester sulfur must be connected (within 10 bonds) to at least one phosphorus (P) atom.
      4. The thioester carbonyl carbon must be bonded to a carbon (other than the sulfur and oxygen)
         that serves as the start of an acyl chain of at least 3 contiguous, aliphatic carbon atoms.
      5. The molecule must contain a purine (adenine) substructure as one indicator of a CoA fragment.
         (Here we try two alternative purine SMARTS patterns to allow for slight variations encountered in 
         the provided structures.)
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): A tuple where the Boolean indicates if the structure qualifies as a fatty acyl-CoA,
                     and the string gives the reason for the decision.
    """
    # Parse the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the thioester fragment: we look for an acyl carbon (C=O) bonded to an S
    thioester_smarts = "C(=O)S"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester (C(=O)S) functional group found"
    
    # Check for a purine (adenine) substructure. We use two alternative SMARTS patterns.
    purine_smarts1 = "n1cnc2ncnc12"  # original pattern (all aromatic)
    purine_smarts2 = "c1nc2c(n1)nc(nc2)N"  # alternative pattern allowing slight differences
    purine_pattern1 = Chem.MolFromSmarts(purine_smarts1)
    purine_pattern2 = Chem.MolFromSmarts(purine_smarts2)
    if not (mol.HasSubstructMatch(purine_pattern1) or mol.HasSubstructMatch(purine_pattern2)):
        return False, "No adenine (purine) substructure found; CoA fragment likely missing"
    
    # Check for presence of phosphorus (P) atoms as part of CoA.
    phosphorus_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not phosphorus_indices:
        return False, "No phosphorus (P) atoms found; CoA moiety likely missing"
    
    reasons = []
    fatty_acyl_found = False

    # Helper: BFS from a start atom to see if a phosphorus atom is reached within threshold bonds.
    def is_connected_to_P(start_idx, threshold=10):
        visited = set([start_idx])
        queue = deque([(start_idx, 0)])
        while queue:
            current, dist = queue.popleft()
            if dist > threshold:
                continue
            if current in phosphorus_indices:
                return True, dist
            for nbr in mol.GetAtomWithIdx(current).GetNeighbors():
                nidx = nbr.GetIdx()
                if nidx not in visited:
                    visited.add(nidx)
                    queue.append((nidx, dist+1))
        return False, None

    # Helper: DFS to count contiguous (aliphatic, non-aromatic) carbon atoms.
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
            # Only count non-aromatic carbons.
            if atom.GetAtomicNum() != 6 or atom.GetIsAromatic():
                continue
            count += 1
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6 and not nbr.GetIsAromatic():
                    stack.append(nbr.GetIdx())
        return count

    # Loop through thioester matches.
    for match in thioester_matches:
        # In "C(=O)S", match[0] → carbonyl carbon, match[1] → sulfur.
        carbonyl_idx = match[0]
        sulfur_idx = match[1]

        # Check connectivity: thioester sulfur must be within 10 bonds of a phosphorus.
        connected, dist = is_connected_to_P(sulfur_idx, threshold=10)
        if not connected:
            reasons.append("Thioester sulfur not connected to a phosphorus within 10 bonds")
            continue

        # From the carbonyl carbon, find a neighbor (other than the sulfur and any oxygen) that is carbon.
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        acyl_start = None
        for nbr in carbonyl_atom.GetNeighbors():
            # Skip the sulfur and oxygens.
            if nbr.GetIdx() == sulfur_idx or nbr.GetAtomicNum() == 8:
                continue
            if nbr.GetAtomicNum() == 6:
                acyl_start = nbr.GetIdx()
                break
        if acyl_start is None:
            reasons.append("No acyl chain found on the carbonyl side")
            continue

        # Count the contiguous aliphatic carbon chain (excluding the carbonyl and sulfur atoms).
        chain_length = count_aliphatic_chain(acyl_start, excluded_idxs={sulfur_idx, carbonyl_idx})
        if chain_length < 3:
            reasons.append(f"Acyl chain too short (chain length {chain_length}); need at least 3 carbon atoms")
            continue

        # If all checks have passed for this thioester fragment, we classify it as a fatty acyl-CoA.
        fatty_acyl_found = True
        break

    if not fatty_acyl_found:
        if reasons:
            return False, "; ".join(reasons)
        else:
            return False, "No valid fatty acyl-CoA structure found"
    
    return True, "Molecule contains a thioester linking a fatty acyl chain with a CoA fragment (purine and P present)"

# Example usage (uncomment to use):
# test_smiles = "S(C(=O)CCC(CCCC(C)C)C)CCNC(=O)CCNC(=O)[C@H](O)C(COP(OP(OC[C@H]1O[C@@H](N2C3=NC=NC(N)=C3N=C2)C(O)[C@H]1OP(O)(O)=O)(O)=O)(O)=O)(C)C"  # Dimethylnonanoyl-CoA
# result, reason = is_fatty_acyl_CoA(test_smiles)
# print(result, reason)