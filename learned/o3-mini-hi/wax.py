"""
Classifies: CHEBI:73702 wax
"""
#!/usr/bin/env python3
"""
Classifies: Wax
Definition: A chemical substance that is an organic compound or mixture of compounds 
           that is composed of long‐chain molecules and is malleable at ambient temperatures.
Approach:
 - Check that the molecule is parsed, has a sufficient molecular weight (≥300 Da) 
   and sufficient rotatable bonds (≥5).
 - Look for ester groups defined by the SMARTS pattern "[CX3](=O)[OX2]".
 - For each ester group found, examine both sides:
     • Acyl (fatty acid) side: from the carbonyl carbon, look at all neighbors (except the linked oxygen) 
       and compute the longest contiguous chain of carbon atoms.
     • Alcohol (fatty alcohol) side: from the oxygen (not linked to the carbonyl), do likewise.
 - The molecule is classified as a wax if there is at least one ester group for which 
   both sides yield a chain of at least 8 carbons.
 - If these criteria are not met, output a reason describing the failure.
  
Note: This is a heuristic and many borderline cases may be misclassified.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_carbon_chain(mol, atom_idx, visited=None):
    """
    Recursively calculates the length of the longest contiguous chain of carbon atoms
    starting from the atom with index atom_idx. Only atoms with atomic number 6 (C) are included.
    
    Args:
        mol (Chem.Mol): The RDKit molecule.
        atom_idx (int): Index of the starting atom.
        visited (set): Set of atom indices that have been visited.
    
    Returns:
        int: The length (number of carbon atoms) in the longest chain from the starting atom.
    """
    if visited is None:
        visited = set()
    visited.add(atom_idx)
    atom = mol.GetAtomWithIdx(atom_idx)
    max_length = 1  # include the starting carbon

    for neighbor in atom.GetNeighbors():
        n_idx = neighbor.GetIdx()
        # continue only if neighbor is carbon and has not been visited yet
        if neighbor.GetAtomicNum() == 6 and n_idx not in visited:
            branch_length = 1 + longest_carbon_chain(mol, n_idx, visited.copy())
            if branch_length > max_length:
                max_length = branch_length
    return max_length

def is_wax(smiles: str):
    """
    Determines if a molecule qualifies as a wax.
    
    The heuristic criteria used are:
      - The molecule must be organic and parsed correctly.
      - It should have a relatively high molecular weight (≥300 Da) and enough rotatable bonds (≥5).
      - It must contain at least one ester group, defined as "[CX3](=O)[OX2]".
      - For at least one ester group, both the acyl (fatty acid) and alcohol (fatty alcohol) parts 
        must have a long contiguous chain (here defined to be ≥8 carbons).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule qualifies as a wax, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Ensure the molecule is organic (has carbon atoms)
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Not organic (contains no carbon atoms)."
    
    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight ({mol_wt:.1f} Da) too low for a typical wax."
    
    # Check number of rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, f"Too few rotatable bonds ({n_rotatable}) for long-chain wax characteristics."
    
    # Define ester SMARTS pattern: carbonyl (C=O) linked to an oxygen
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester functional groups found (required for wax classification)."
    
    # Iterate over each ester match and evaluate chain lengths on both sides.
    for match in ester_matches:
        # match[0] is the carbonyl carbon, match[1] is the oxygen in the ester bond.
        carbonyl_idx = match[0]
        oxygen_idx = match[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        oxygen_atom = mol.GetAtomWithIdx(oxygen_idx)
        
        # For the acyl (fatty acid) side: check all neighbors of the carbonyl carbon except the ester oxygen.
        acid_chain_length = 0
        for neighbor in carbonyl_atom.GetNeighbors():
            if neighbor.GetIdx() == oxygen_idx:
                continue
            if neighbor.GetAtomicNum() == 6:
                candidate_length = longest_carbon_chain(mol, neighbor.GetIdx())
                if candidate_length > acid_chain_length:
                    acid_chain_length = candidate_length
        
        # For the alcohol (fatty alcohol) side: check all neighbors of the ester oxygen except the carbonyl carbon.
        alcohol_chain_length = 0
        for neighbor in oxygen_atom.GetNeighbors():
            if neighbor.GetIdx() == carbonyl_idx:
                continue
            if neighbor.GetAtomicNum() == 6:
                candidate_length = longest_carbon_chain(mol, neighbor.GetIdx())
                if candidate_length > alcohol_chain_length:
                    alcohol_chain_length = candidate_length
        
        if acid_chain_length >= 8 and alcohol_chain_length >= 8:
            return True, (f"Found an ester group with acyl chain length {acid_chain_length} "
                          f"and alcohol chain length {alcohol_chain_length}, consistent with a wax.")
    
    return False, "No ester group found with two sufficiently long (≥8 carbons) chains."

# Example usage:
if __name__ == "__main__":
    # Example: 2-palmitoyloxypalmityl palmitate (one of the wax examples)
    test_smiles = "CCCCCCCCCCCCCCCC(=O)OCC(CCCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCCCC"
    classification, reason = is_wax(test_smiles)
    print(f"Classification: {classification}\nReason: {reason}")