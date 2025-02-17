"""
Classifies: CHEBI:73702 wax
"""
#!/usr/bin/env python3
"""
Classifies: Wax
Definition: A chemical substance that is an organic compound or mixture of compounds 
           that is composed of long-chain molecules and is malleable at ambient temperatures.
Approach:
 - First, check that the molecule is organic (has carbon atoms), has sufficient molecular weight (>=300 Da)
   and enough flexibility (>=5 rotatable bonds).
 - Search for ester groups using the SMARTS pattern "[CX3](=O)[OX2]".
 - For each ester group, look at the two sides of the ester:
       • For the acyl side (from the carbonyl carbon), examine all neighbors except the ester oxygen.
       • For the alcohol side (from the ester oxygen), examine all neighbors except the carbonyl carbon.
   Use a backtracking depth‐first search (DFS) to calculate the longest contiguous chain of carbon atoms.
 - If at least one ester group is found where both the acyl and alcohol chains have a longest contiguous chain 
   length of at least 8 carbons, we classify the molecule as a wax.
Note: This heuristic is approximate.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_chain_from(mol, current_idx, exclude: set, visited: set) -> int:
    """
    Backtracking DFS to compute the longest contiguous chain of carbon atoms
    starting from the atom given by current_idx.
    Args:
        mol (Chem.Mol): RDKit molecule.
        current_idx (int): Starting atom index (should be carbon).
        exclude (set): Atom indices that should not be visited.
        visited (set): Atom indices already visited along this path.
    Returns:
        int: The length (number of carbons) of the longest chain from the starting atom.
    """
    visited.add(current_idx)
    max_length = 1  # Count the current carbon
    atom = mol.GetAtomWithIdx(current_idx)
    for nbr in atom.GetNeighbors():
        nbr_idx = nbr.GetIdx()
        if nbr_idx in visited or nbr_idx in exclude:
            continue
        # Only traverse into carbons (atomic number 6)
        if nbr.GetAtomicNum() == 6:
            length = 1 + longest_chain_from(mol, nbr_idx, exclude, visited)
            if length > max_length:
                max_length = length
    visited.remove(current_idx)
    return max_length

def is_wax(smiles: str):
    """
    Determines if a molecule qualifies as a wax.
    
    Heuristic criteria:
      - Must be an organic molecule (contain carbon atoms)
      - Must have a molecular weight ≥ 300 Da and at least 5 rotatable bonds.
      - Must contain at least one ester group "[CX3](=O)[OX2]" where, after “cutting” the ester bond,
        both the acyl (fatty acid) part and the alcohol (fatty alcohol) part have a longest 
        contiguous carbon chain length of at least 8.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule qualifies as wax, False otherwise.
        str: Detailed reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    # Sanitize the molecule to ensure proper valence and connectivity.
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization failed: {e}"
    
    # Check that the molecule is organic (must have carbon atoms)
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Not organic (no carbon atoms found)."
    
    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight ({mol_wt:.1f} Da) too low for a typical wax."
    
    # Check rotatable bonds count
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, f"Too few rotatable bonds ({n_rotatable}) for long-chain wax characteristics."
    
    # Define and find ester groups.
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester functional groups found (required for wax classification)."
    
    # Loop over each ester match
    for match in ester_matches:
        # In the matched ester, match[0] is the carbonyl carbon and match[1] is the ester oxygen.
        carbonyl_idx = match[0]
        oxygen_idx = match[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        oxygen_atom = mol.GetAtomWithIdx(oxygen_idx)
        
        # For acyl (fatty acid) side: examine neighbors of carbonyl carbon excluding the ester oxygen.
        acyl_chain_length = 0
        acyl_valid = False
        for nbr in carbonyl_atom.GetNeighbors():
            if nbr.GetIdx() == oxygen_idx:
                continue
            if nbr.GetAtomicNum() == 6:  # must be carbon
                # Compute the longest contiguous chain.
                chain_len = longest_chain_from(mol, nbr.GetIdx(), exclude={oxygen_idx}, visited=set())
                # Optionally, add 1 to include the neighbor itself (already done in our DFS)
                if chain_len >= 8:
                    acyl_chain_length = max(acyl_chain_length, chain_len)
                    acyl_valid = True
        # For alcoholic (fatty alcohol) side: examine neighbors of oxygen excluding the carbonyl carbon.
        alcohol_chain_length = 0
        alcohol_valid = False
        for nbr in oxygen_atom.GetNeighbors():
            if nbr.GetIdx() == carbonyl_idx:
                continue
            if nbr.GetAtomicNum() == 6:
                chain_len = longest_chain_from(mol, nbr.GetIdx(), exclude={carbonyl_idx}, visited=set())
                if chain_len >= 8:
                    alcohol_chain_length = max(alcohol_chain_length, chain_len)
                    alcohol_valid = True
        
        # Check if both sides have a sufficiently long chain.
        if acyl_valid and alcohol_valid:
            return True, (f"Found an ester group with acyl chain length {acyl_chain_length} and "
                          f"alcohol chain length {alcohol_chain_length}, consistent with wax compounds.")
    
    return False, "No ester group found with two sufficiently long (≥8 carbons) chains."

# Example usage:
if __name__ == "__main__":
    # Test with one of the provided SMILES: 2-palmitoyloxypalmityl palmitate
    test_smiles = "CCCCCCCCCCCCCCCC(=O)OCC(CCCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCCCC"
    classification, reason = is_wax(test_smiles)
    print(f"Classification: {classification}\nReason: {reason}")