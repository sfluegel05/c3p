"""
Classifies: CHEBI:90546 medium-chain fatty acyl-CoA(4-)
"""
#!/usr/bin/env python
"""
Classifies: medium-chain fatty acyl-CoA(4-)
Definition: An acyl-CoA oxoanion that results from deprotonation of the phosphate and diphosphate groups 
of any medium-chain fatty acyl-CoA; major species at pH 7.3.
The criteria here are:
  1. SMILES must be valid.
  2. The molecule must contain a thioester bond (a carbonyl bound to sulfur).
  3. It must contain two characteristic substructures of Coenzyme A:
       - an adenine (nucleobase) ring,
       - and a pantetheine-like fragment (SCCNC(=O)CCNC(=O)) present in most acyl‐CoAs.
  4. The overall formal charge must be -4.
  5. The acyl chain attached at the thioester must be “medium‐chain” in that,
     when counting only “linear” (acyclic) carbon atoms (starting at the carbonyl carbon),
     the chain length (including the carbonyl carbon) is between 6 and 12.
     
If any of these criteria fail, we return False and an explanation.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_medium_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a medium-chain fatty acyl-CoA(4-), False otherwise.
        str: Explanation for the decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Look for the thioester bond:
    # We use a SMARTS that identifies a carbonyl carbon bound to a sulfur.
    thioester_smarts = "[CX3](=O)S"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester (acyl-CoA) bond found"
    
    # 2. Check for CoA-like fragments.
    # (a) The nucleobase part (adenine-type) and 
    # (b) a portion of the pantetheine moiety.
    adenine_pattern = Chem.MolFromSmarts("n1cnc2c(ncnc12)")
    pantetheine_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "No adenine (nucleobase) fragment found"
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "No pantetheine fragment (SCCNC(=O)CCNC(=O)) found"
    
    # 3. Check overall formal charge.
    net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if net_charge != -4:
        return False, f"Net charge is {net_charge}, but expected -4 for this class"
    
    # 4. Extract the acyl chain.
    # We pick one thioester match (usually there is only one).
    # In the match tuple, the first atom (index 0) corresponds to the carbonyl carbon.
    carbonyl_idx = thioester_matches[0][0]
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    
    # Now, among the neighbors of the carbonyl carbon, we need to pick the one that is part of the acyl chain.
    # (We discard the sulfur neighbor.)
    acyl_neighbors = [nbr for nbr in carbonyl_atom.GetNeighbors() if nbr.GetSymbol() == "C"]
    if not acyl_neighbors:
        return False, "No acyl chain attached to the thioester bond"
    
    # Define a DFS routine to compute the maximum number of contiguous acyclic (non-ring) carbons,
    # starting from a given carbon atom.
    def dfs_chain_length(atom, visited):
        # Only traverse through carbon atoms that are not in any ring.
        current_idx = atom.GetIdx()
        visited.add(current_idx)
        max_len = 1  # counting current carbon atom
        for nbr in atom.GetNeighbors():
            if nbr.GetSymbol() != "C":
                continue
            # Only follow carbon atoms that are not part of any ring.
            if nbr.IsInRing():
                continue
            nbr_idx = nbr.GetIdx()
            if nbr_idx in visited:
                continue
            length = 1 + dfs_chain_length(nbr, visited.copy())
            if length > max_len:
                max_len = length
        return max_len

    # For each acyl neighbor, calculate the chain length.
    # We include the carbonyl carbon as part of the chain.
    candidate_chain_length = 0
    for acyl_atom in acyl_neighbors:
        # Only consider neighbor if it is not in a ring (to enforce a linear chain).
        if acyl_atom.IsInRing():
            continue
        chain_len = 1 + dfs_chain_length(acyl_atom, set([carbonyl_idx]))
        if chain_len > candidate_chain_length:
            candidate_chain_length = chain_len

    if candidate_chain_length == 0:
        return False, "Could not extract an acyl chain from the thioester bond"

    # 5. Check that the acyl chain length is within the medium-chain range (6 to 12 carbons, including the carbonyl).
    if candidate_chain_length < 6 or candidate_chain_length > 12:
        return False, f"Acyl chain length is {candidate_chain_length} carbons, not in medium-chain range (6-12)"
    
    return True, (f"Matches acyl-CoA pattern with medium chain length of {candidate_chain_length} carbons "
                  "and net charge -4")

# Example usage:
if __name__ == "__main__":
    # Test with one example SMILES string for cis-dec-3-enoyl-CoA(4-)
    test_smiles = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(=O)C/C=C\\CCCCCC)=O)=O)O)C)(=O)[O-])(=O)[O-])[C@H]([C@H]1O)OP([O-])([O-])=O"
    result, reason = is_medium_chain_fatty_acyl_CoA_4__(test_smiles)
    print(result, reason)