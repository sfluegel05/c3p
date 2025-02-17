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
    Our criteria:
      1. The SMILES must be valid.
      2. Contains a thioester group ([C;X3](=O)[S]).
      3. The thioester sulfur must be within 15 bonds of at least one phosphorus (P) atom.
      4. The carbonyl carbon of the thioester must lead to an aliphatic acyl chain (≥3 contiguous non-aromatic carbons).
      5. The molecule must contain at least one purine (adenine) substructure, detected via several SMARTS variants.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule qualifies as a fatty acyl-CoA, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Look for the thioester functional group.
    # Pattern: a carbon (with three connections) double-bonded to an oxygen and singly bound to a sulfur.
    thioester_pattern = Chem.MolFromSmarts("[C;X3](=O)[S]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester (C(=O)S) functional group found"
    
    # 2. Look for a purine (adenine) substructure.
    # We include multiple SMARTS: some that use aromatic (lowercase) notation plus an explicit uppercase variant.
    purine_smarts = [
        "n1cnc2ncnc12",          # common aromatic adenine pattern
        "N1C=NC2=C1N=CN2",        # explicit adenine pattern (uppercase)
        "c1nc2c(n1)nc(nc2)N",      # alternative aromatic pattern
        "n1c([nH])nc2ncnc12"       # additional variant
    ]
    purine_patterns = [Chem.MolFromSmarts(s) for s in purine_smarts if Chem.MolFromSmarts(s) is not None]
    if not any(mol.HasSubstructMatch(p) for p in purine_patterns):
        return False, "No adenine (purine) substructure found; CoA fragment likely missing"
    
    # 3. Check that there is at least one phosphorus (P) atom.
    phosphorus_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not phosphorus_indices:
        return False, "No phosphorus (P) atoms found; CoA moiety likely missing"
    
    reasons = []
    fatty_acyl_found = False
    
    # Helper: BFS from a starting atom index to determine if a phosphorus atom is connected within threshold bonds.
    def is_connected_to_P(start_idx, threshold=15):
        visited = set([start_idx])
        queue = deque([(start_idx, 0)])
        while queue:
            current, dist = queue.popleft()
            if current in phosphorus_indices:
                return True, dist
            if dist >= threshold:
                continue
            for nbr in mol.GetAtomWithIdx(current).GetNeighbors():
                nidx = nbr.GetIdx()
                if nidx not in visited:
                    visited.add(nidx)
                    queue.append((nidx, dist + 1))
        return False, None
    
    # Helper: Count contiguous, non-aromatic carbon atoms starting from a given index (exclude certain indices).
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

    # 4. Process each thioester match.
    for match in thioester_matches:
        # In the pattern "[C](=O)[S]", match[0] is the carbonyl carbon and match[1] is the sulfur.
        carbonyl_idx = match[0]
        sulfur_idx = match[1]
        
        # 4a. Check that thioester sulfur is connected within 15 bonds to a phosphorus atom.
        connected, dist = is_connected_to_P(sulfur_idx, threshold=15)
        if not connected:
            reasons.append("Thioester sulfur not connected to a phosphorus within 15 bonds")
            continue
        
        # 4b. From the carbonyl carbon, find a neighbor (other than sulfur or oxygen) that is a carbon atom.
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        acyl_start = None
        for nbr in carbonyl_atom.GetNeighbors():
            if nbr.GetIdx() == sulfur_idx or nbr.GetAtomicNum() == 8:
                continue
            if nbr.GetAtomicNum() == 6:
                acyl_start = nbr.GetIdx()
                break
        if acyl_start is None:
            reasons.append("No acyl chain found on the carbonyl side")
            continue
        
        # 4c. Count the contiguous aliphatic carbon chain (excluding the sulfur and the carbonyl carbon).
        chain_length = count_aliphatic_chain(acyl_start, excluded_idxs={sulfur_idx, carbonyl_idx})
        if chain_length < 3:
            reasons.append(f"Acyl chain too short (chain length {chain_length}); need at least 3 carbon atoms")
            continue
        
        # This thioester fragment passes all tests.
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