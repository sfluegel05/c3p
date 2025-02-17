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
    
    # Check for the thioester bond with a carbonyl bound to sulfur.
    # Pattern: a carbon (with a double bond to oxygen) bonded to a sulfur.
    thioester_pattern = Chem.MolFromSmarts("[#6](=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester (acyl-CoA) bond found"
    
    # Look for an adenine-like fragment to mark the CoA portion.
    coa_pattern = Chem.MolFromSmarts("n1cnc2c(ncnc12)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA (adenine/nucleoside) fragment found"
    
    # Check the overall formal charge. For this class we expect a net charge of -4.
    net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if net_charge != -4:
        return False, f"Net charge is {net_charge}, but expected -4 for this class"
    
    # Now, isolate the fatty acyl chain.
    # We use one of the thioester matches (usually there is only one).
    # The thioester match is a tuple of atom indices; we take the first index
    # as the carbonyl carbon and find the connected carbon (not the sulfur) that forms the acyl chain.
    chain_found = False
    chain_length = 0
    for match in thioester_matches:
        carbonyl_idx = match[0]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # Identify the neighbor atoms of the carbonyl carbon that are carbons.
        acyl_neighbors = [nbr for nbr in carbonyl_atom.GetNeighbors() if nbr.GetSymbol() == "C"]
        if not acyl_neighbors:
            continue  # try next match if exists

        # Define a DFS routine that tracks visited atoms to account for cycles.
        def dfs(atom, visited):
            """Return the longest chain length (number of carbons) starting from atom,
               while avoiding cycles.
               visited is a set of atom indices already being traversed."""
            current_idx = atom.GetIdx()
            visited.add(current_idx)
            max_len = 1  # Count current carbon atom.
            for nbr in atom.GetNeighbors():
                if nbr.GetSymbol() != "C":
                    continue
                nbr_idx = nbr.GetIdx()
                if nbr_idx in visited:
                    continue
                length = 1 + dfs(nbr, visited.copy())  # copy visited for independent branch traversal.
                if length > max_len:
                    max_len = length
            return max_len
        
        # For each acyl neighbor, calculate chain length and pick the longest.
        candidate_max = 0
        for acyl_atom in acyl_neighbors:
            # Include the carbonyl carbon.
            path_len = 1 + dfs(acyl_atom, set([carbonyl_idx]))
            if path_len > candidate_max:
                candidate_max = path_len
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
    # Test with one example SMILES string
    test_smiles = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(=O)C/C=C\\CCCCCC)=O)=O)O)C)(=O)[O-])(=O)[O-])[C@H]([C@H]1O)OP([O-])([O-])=O"
    result, reason = is_medium_chain_fatty_acyl_CoA_4__(test_smiles)
    print(result, reason)