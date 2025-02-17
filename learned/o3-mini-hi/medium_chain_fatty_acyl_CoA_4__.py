"""
Classifies: CHEBI:90546 medium-chain fatty acyl-CoA(4-)
"""
#!/usr/bin/env python
"""
Classifies: medium-chain fatty acyl-CoA(4-)
Definition: An acyl-CoA oxoanion that results from deprotonation of the phosphate and diphosphate groups 
of any medium-chain fatty acyl-CoA; major species at pH 7.3.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_medium_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA(4-)
    based on its SMILES string.
    
    We require that:
      1. The SMILES is valid.
      2. It contains the thioester linkage (a carbonyl bound to sulfur) that marks the acyl-CoA bond.
      3. It contains a CoA-like moiety (here identified by an adenine/nucleoside substructure).
      4. The net formal charge of the molecule is -4.
      5. The fatty acyl (chain) part, taken from the thioester carbonyl atom, 
         has a linear chain between 6 and 12 carbon atoms (inclusive).
         
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as a medium-chain fatty acyl-CoA(4-), False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the thioester bond with a carbonyl: pattern [#6](=O)S
    thioester_pattern = Chem.MolFromSmarts("[#6](=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester (acyl-CoA) bond found"
    
    # For CoA moiety we look for an adenine fragment.
    # Many CoA compounds include a nucleoside portion containing adenine.
    coa_pattern = Chem.MolFromSmarts("n1cnc2c(ncnc12)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA (adenine/nucleoside) fragment found"
    
    # Check the overall formal charge. For this class we expect a net charge of -4.
    net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if net_charge != -4:
        return False, f"Net charge is {net_charge}, but expected -4 for this class"
    
    # Now, try to isolate the fatty acyl chain length.
    # We look at one of the thioester matches. (In many acyl-CoAs there is only one thioester.)
    # thioester Match: tuple of atom indices; the first index is the carbonyl carbon.
    # From the carbonyl carbon, its neighbors include oxygen (of the carbonyl) and S.
    # The acyl chain is the chain attached to the carbonyl (not the S side).
    chain_found = False
    chain_length = 0
    for match in thioester_matches:
        carbonyl_idx = match[0]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # find neighbor that is a carbon (this is the acyl chain substituent)
        acyl_neighbors = [nbr for nbr in carbonyl_atom.GetNeighbors() 
                          if nbr.GetSymbol() == "C"]
        if not acyl_neighbors:
            continue  # try next match if exists
        # choose the first candidate (or the longest path among candidates)
        # We will count the chain length using a DFS on carbon atoms.
        # Starting count includes the carbonyl carbon as part of the chain.
        def dfs(atom, came_from_idx):
            """Return the longest path length (number of carbons) starting from atom.
               Only follow bonds between carbon atoms."""
            max_len = 1  # Count current atom
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() == came_from_idx:
                    continue
                if nbr.GetSymbol() != "C":
                    continue
                # Continue DFS on neighbor.
                length = 1 + dfs(nbr, atom.GetIdx())
                if length > max_len:
                    max_len = length
            return max_len
        
        # Get the maximum chain length by considering each acyl neighbor.
        candidate_max = 0
        for acyl_atom in acyl_neighbors:
            path_len = 1 + dfs(acyl_atom, carbonyl_idx)  # plus the carbonyl carbon
            if path_len > candidate_max:
                candidate_max = path_len
        # If a chain was found, record it and break out.
        if candidate_max:
            chain_found = True
            chain_length = candidate_max
            break
        
    if not chain_found:
        return False, "No acyl chain could be extracted from the thioester bond"
    
    # Typically, medium-chain fatty acyl groups have between 6 and 12 carbons (including the carbonyl carbon).
    if chain_length < 6 or chain_length > 12:
        return False, f"Acyl chain length is {chain_length} carbons, not in medium-chain range (6-12)"
    
    return True, f"Matches acyl-CoA pattern with medium chain length of {chain_length} carbons and net charge -4"

# Example usage:
if __name__ == "__main__":
    # Test with one example SMILES (e.g., cis-dec-3-enoyl-CoA(4-))
    test_smiles = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(=O)C/C=C\\CCCCCC)=O)=O)O)C)(=O)[O-])(=O)[O-])[C@H]([C@H]1O)OP([O-])([O-])=O"
    result, reason = is_medium_chain_fatty_acyl_CoA_4__(test_smiles)
    print(result, reason)