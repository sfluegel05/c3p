"""
Classifies: CHEBI:37554 fatty acyl-CoA
"""
"""
Classifies: fatty acyl-CoA
Definition: An acyl-CoA results from the condensation of the thiol group of coenzyme A with the carboxy group of a fatty acid.
The molecule must contain a thioester group (C(=O)S) linking an acyl chain on the carbonyl side to a CoA fragment.
We approximate the CoA fragment by requiring (a) at least one phosphorus (P) atom and (b) a purine (adenine) substructure.
"""

from rdkit import Chem
from collections import deque

def is_fatty_acyl_CoA(smiles: str):
    """
    Determines if the given SMILES string corresponds to a fatty acyl-CoA.
    Our criteria are:
      1. The SMILES is valid.
      2. Must contain a thioester group (C(=O)S).
      3. The thioester sulfur must be connected (within 15 bonds) to at least one phosphorus (P) atom.
      4. The thioester carbonyl carbon must be bonded to a carbon (other than sulfur or oxygen)
         that initiates an aliphatic chain of at least 3 contiguous, non-aromatic carbon atoms.
      5. The molecule must contain a purine (adenine) substructure as a marker of a CoA fragment.
         Here we try three alternative purine SMARTS patterns.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): A tuple with a boolean (True if the structure qualifies as a fatty acyl-CoA)
                     and a string with the reason for the classification.
    """
    # Parse the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for the thioester fragment.
    # Slightly modified to catch possible stereochemistry details.
    thioester_smarts = "[C;X3](=O)[S]"  # Acyl carbonyl (with three connections) attached to an S.
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester (C(=O)S) functional group found"
    
    # 2. Check for a purine (adenine) substructure.
    purine_smarts1 = "n1cnc2ncnc12"  # common adenine pattern
    purine_smarts2 = "c1nc2c(n1)nc(nc2)N"  # alternative adenine pattern
    purine_smarts3 = "n1c([nH])nc2ncnc12"  # additional variant
    purine_patterns = [Chem.MolFromSmarts(s) for s in [purine_smarts1, purine_smarts2, purine_smarts3]]
    if not any(mol.HasSubstructMatch(p) for p in purine_patterns):
        return False, "No adenine (purine) substructure found; CoA fragment likely missing"
    
    # 3. Check for at least one phosphorus (P) atom in the structure.
    phosphorus_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not phosphorus_indices:
        return False, "No phosphorus (P) atoms found; CoA moiety likely missing"
    
    reasons = []
    fatty_acyl_found = False

    # Helper: Breadth-first search (BFS) from a start atom index to see if a phosphorus atom is reached
    # within a threshold number of bonds.
    def is_connected_to_P(start_idx, threshold=15):
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
                    queue.append((nidx, dist + 1))
        return False, None

    # Helper: Count contiguous aliphatic (non-aromatic) carbon atoms starting from a given carbon.
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
            # Count only non-aromatic carbons.
            if atom.GetAtomicNum() != 6 or atom.GetIsAromatic():
                continue
            count += 1
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6 and not nbr.GetIsAromatic():
                    stack.append(nbr.GetIdx())
        return count

    # 4. Loop through all thioester matches.
    for match in thioester_matches:
        # For the pattern "[C](=O)[S]", match[0] → carbonyl carbon, match[1] → sulfur.
        carbonyl_idx = match[0]
        sulfur_idx = match[1]

        # 4a. Check connectivity: thioester sulfur must be within threshold bonds of a phosphorus atom.
        connected, dist = is_connected_to_P(sulfur_idx, threshold=15)
        if not connected:
            reasons.append("Thioester sulfur not connected to a phosphorus within 15 bonds")
            continue  # Try the next thioester match if any

        # 4b. From the carbonyl carbon, look for a neighbor (other than the sulfur or oxygen) that is carbon:
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        acyl_start = None
        for nbr in carbonyl_atom.GetNeighbors():
            # Skip the sulfur and any oxygens.
            if nbr.GetIdx() == sulfur_idx or nbr.GetAtomicNum() == 8:
                continue
            if nbr.GetAtomicNum() == 6:
                acyl_start = nbr.GetIdx()
                break
        if acyl_start is None:
            reasons.append("No acyl chain found on the carbonyl side")
            continue

        # 4c. Count the contiguous aliphatic carbon chain (exclude the carbonyl and sulfur atoms).
        chain_length = count_aliphatic_chain(acyl_start, excluded_idxs={sulfur_idx, carbonyl_idx})
        if chain_length < 3:
            reasons.append(f"Acyl chain too short (chain length {chain_length}); need at least 3 carbon atoms")
            continue

        # If all checks have been passed for this thioester fragment, we classify it as fatty acyl-CoA.
        fatty_acyl_found = True
        break

    if not fatty_acyl_found:
        if reasons:
            return False, "; ".join(reasons)
        else:
            return False, "No valid fatty acyl-CoA structure found"
    
    return True, "Molecule contains a thioester linking a fatty acyl chain with a CoA fragment (purine and P present)"

# Example usage:
# test_smiles = "S(C(=O)CCC(CCCC(C)C)C)CCNC(=O)CCNC(=O)[C@H](O)C(COP(OP(OC[C@H]1O[C@@H](N2C3=NC=NC(N)=C3N=C2)C(O)[C@H]1OP(O)(O)=O)(O)=O)(O)=O)(C)C"  # Dimethylnonanoyl-CoA
# result, reason = is_fatty_acyl_CoA(test_smiles)
# print(result, reason)