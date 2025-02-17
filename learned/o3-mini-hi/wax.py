"""
Classifies: CHEBI:73702 wax
"""
#!/usr/bin/env python3
"""
Classifies: Wax
Definition: A chemical substance that is an organic compound or mixture of compounds 
           that is composed of long‐chain molecules and is malleable at ambient temperatures.
Approach:
 - Parse the molecule and check that it is organic (has carbons), has a high enough 
   molecular weight (≥300 Da) and a sufficient number of rotatable bonds (≥5).
 - Find ester groups using the SMARTS "[CX3](=O)[OX2]".
 - For each ester group, separate the two sides:
       • For the acyl (fatty acid) part: from the carbonyl carbon, 
         ignore the bond to the ester oxygen and perform a DFS to count the longest contiguous carbon chain.
       • For the alcohol (fatty alcohol) part: from the ester oxygen, 
         ignore the bond back to the carbonyl carbon and perform a DFS to count the longest contiguous carbon chain.
   The DFS is implemented with an exclude set to avoid crossing back over the ester bond.
 - Return True if any ester group has both chains ≥8 carbons.
Note: This heuristic may be imperfect but should catch the majority of expected wax compounds.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_carbon_chain_from(mol, start_idx, exclude: set = None, visited: set = None) -> int:
    """
    Recursively count the longest chain of connected carbon atoms starting at start_idx,
    while not crossing into any atom index in the exclude set.
    
    Args:
        mol (Chem.Mol): RDKit molecule.
        start_idx (int): Starting atom index (should be carbon).
        exclude (set): set of atom indices to not traverse.
        visited (set): set of already visited atom indices (for this DFS branch).
        
    Returns:
        int: Length (number of carbons) of the longest contiguous chain from start_idx.
    """
    if exclude is None:
        exclude = set()
    if visited is None:
        visited = set()
    visited.add(start_idx)
    max_length = 1  # count the current carbon

    atom = mol.GetAtomWithIdx(start_idx)
    for nbr in atom.GetNeighbors():
        nbr_idx = nbr.GetIdx()
        # Only continue if neighbor is carbon, not in visited, and not excluded.
        if nbr_idx in visited or nbr_idx in exclude:
            continue
        if nbr.GetAtomicNum() == 6:
            branch_length = 1 + longest_carbon_chain_from(mol, nbr_idx, exclude, visited.copy())
            if branch_length > max_length:
                max_length = branch_length
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
        str: Reason for classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Check organic: must contain carbon atoms.
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Not organic (contains no carbon atoms)."
    
    # Check molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight ({mol_wt:.1f} Da) too low for a typical wax."
    
    # Check rotatable bonds.
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, f"Too few rotatable bonds ({n_rotatable}) for long-chain wax characteristics."
    
    # Define ester SMARTS pattern.
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester functional groups found (required for wax classification)."
    
    # Loop over each ester group.
    for match in ester_matches:
        # match[0]: carbonyl carbon; match[1]: ester oxygen.
        carbonyl_idx = match[0]
        oxygen_idx = match[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        oxygen_atom = mol.GetAtomWithIdx(oxygen_idx)
        
        # For the acyl (fatty acid) side: 
        # select neighbors of the carbonyl carbon except the oxygen of the ester.
        acyl_chain_length = 0
        for neighbor in carbonyl_atom.GetNeighbors():
            if neighbor.GetIdx() == oxygen_idx:
                continue  # exclude the ester oxygen
            if neighbor.GetAtomicNum() == 6:
                # Exclude the carbonyl carbon itself from further traversal to avoid jumping back.
                chain_length = longest_carbon_chain_from(mol, neighbor.GetIdx(), exclude={oxygen_idx})
                if chain_length > acyl_chain_length:
                    acyl_chain_length = chain_length
        
        # For the alcoholic (fatty alcohol) side:
        # select neighbors of the oxygen atom except the carbonyl carbon.
        alcohol_chain_length = 0
        for neighbor in oxygen_atom.GetNeighbors():
            if neighbor.GetIdx() == carbonyl_idx:
                continue  # exclude the carbonyl carbon
            if neighbor.GetAtomicNum() == 6:
                chain_length = longest_carbon_chain_from(mol, neighbor.GetIdx(), exclude={carbonyl_idx})
                if chain_length > alcohol_chain_length:
                    alcohol_chain_length = chain_length

        # Check if both chains are sufficiently long (≥8 carbons).
        if acyl_chain_length >= 8 and alcohol_chain_length >= 8:
            return True, (f"Found an ester group with acyl chain length {acyl_chain_length} "
                          f"and alcohol chain length {alcohol_chain_length}, consistent with a wax.")
    
    return False, "No ester group found with two sufficiently long (≥8 carbons) chains."

# Example usage
if __name__ == "__main__":
    test_smiles = "CCCCCCCCCCCCCCCC(=O)OCC(CCCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCCCC"  # 2-palmitoyloxypalmityl palmitate
    classification, reason = is_wax(test_smiles)
    print(f"Classification: {classification}\nReason: {reason}")