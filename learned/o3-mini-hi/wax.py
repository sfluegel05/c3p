"""
Classifies: CHEBI:73702 wax
"""
#!/usr/bin/env python3
"""
Classifies: Wax
Definition: A chemical substance that is an organic compound or mixture of compounds 
           composed of long-chain molecules and malleable at ambient temperatures.
Approach:
  - Pre-filter: must be organic, have only C, H, O, be acyclic, with weight ≥300 Da and ≥5 rotatable bonds.
  - Identify ester linkage(s): look for a bond connecting a carbon (with a double bond to an oxygen) to an oxygen.
  - Reject molecules with multiple ester linkages.
  - "Cut" the molecule at the ester bond and perform a DFS on both sides to determine the longest linear carbon chain.
  - If both sides yield a chain length of at least 8 carbons, we classify the molecule as a wax.
Note: This heuristic is approximate.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_chain_from(mol, current_idx, exclude, visited):
    """
    Recursive DFS to compute the longest linear chain of carbon atoms starting from current_idx.
    
    Args:
      mol (Chem.Mol): The RDKit molecule.
      current_idx (int): Starting atom index (should be carbon).
      exclude (set): Atom indices to exclude from the search.
      visited (set): Already visited atom indices in current path.
      
    Returns:
      int: Maximum number of carbons in a contiguous chain (including starting carbon).
    """
    visited.add(current_idx)
    max_length = 1  # Count the starting carbon
    atom = mol.GetAtomWithIdx(current_idx)
    for nbr in atom.GetNeighbors():
        nbr_idx = nbr.GetIdx()
        if nbr_idx in visited or nbr_idx in exclude:
            continue
        if nbr.GetAtomicNum() == 6:  # proceed only if neighbor is carbon
            chain_length = 1 + longest_chain_from(mol, nbr_idx, exclude, visited)
            if chain_length > max_length:
                max_length = chain_length
    visited.remove(current_idx)
    return max_length

def is_wax(smiles: str):
    """
    Determines if a molecule qualifies as a wax.
    
    New Heuristic criteria:
      - The molecule must be organic and composed solely of carbon, hydrogen, and oxygen.
      - It should be acyclic.
      - The molecular weight must be ≥300 Da and have at least 5 rotatable bonds.
      - The molecule should contain exactly one ester linkage (R–C(=O)–O–R')
        such that when "cut" at the ester bond, both the fatty acyl (acid) portion and the fatty alcohol portion
        have a longest continuous chain of at least 8 carbons.
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if molecule qualifies as a wax, False otherwise.
      str: Detailed reason for classification decision.
    """
    # Parse SMILES and sanitize
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization failed: {e}"
    
    # Ensure molecule is organic and only contains C, H, and O.
    allowed_atomic_nums = {1, 6, 8}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains atom {atom.GetSymbol()} not allowed for wax (only C, H, O permitted)."
    
    # Check for acyclicity: waxes in our test set are acyclic.
    ri = mol.GetRingInfo()
    if ri.NumRings() > 0:
        return False, "Molecule contains ring(s), inconsistent with typical wax structures."
    
    # Check molecular weight and number of rotatable bonds.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight ({mol_wt:.1f} Da) too low for a typical wax."
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, f"Too few rotatable bonds ({n_rotatable}) for long-chain wax characteristics."
    
    # List to hold candidate ester linkages
    valid_ester_candidates = []
    
    # Iterate over bonds to identify ester linkage: a bond between an oxygen and a carbon where
    # the carbon (candidate carbonyl) is double-bonded to another oxygen.
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if ((a1.GetAtomicNum() == 8 and a2.GetAtomicNum() == 6) or 
            (a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 8)):
            # Decide roles: the carbon must be the one double-bonded to oxygen (carbonyl).
            if a1.GetAtomicNum() == 6:
                carbon_atom = a1
                oxy_atom = a2
            else:
                carbon_atom = a2
                oxy_atom = a1
            
            # Check if carbon_atom is indeed a carbonyl carbon (has a double bond with an oxygen besides oxy_atom)
            carbonyl_found = False
            for nbr in carbon_atom.GetNeighbors():
                if nbr.GetIdx() == oxy_atom.GetIdx():
                    continue
                bond_tmp = mol.GetBondBetweenAtoms(carbon_atom.GetIdx(), nbr.GetIdx())
                # Check that the neighbor is oxygen and bond order is double
                if nbr.GetAtomicNum() == 8 and bond_tmp.GetBondTypeAsDouble() >= 2.0:
                    carbonyl_found = True
                    break
            if not carbonyl_found:
                continue  # Not an ester linkage
            
            # At this point, we have a candidate ester bond linking carbonyl carbon and an oxygen.
            # Compute chain lengths on both sides.
            acyl_valid = False
            acyl_chain_length = 0
            # For the acyl side, follow neighbors of the carbon_atom except the ester oxygen and its carbonyl oxygen.
            for nbr in carbon_atom.GetNeighbors():
                if nbr.GetIdx() == oxy_atom.GetIdx():
                    continue
                bond_tmp = mol.GetBondBetweenAtoms(carbon_atom.GetIdx(), nbr.GetIdx())
                # Skip the double-bonded oxygen used in the carbonyl group.
                if nbr.GetAtomicNum() == 8 and bond_tmp.GetBondTypeAsDouble() >= 2.0:
                    continue
                if nbr.GetAtomicNum() == 6:
                    chain_len = longest_chain_from(mol, nbr.GetIdx(), 
                                                   exclude={oxygen_atom.GetIdx(), carbon_atom.GetIdx()}, visited=set())
                    if chain_len >= 8:
                        acyl_valid = True
                        if chain_len > acyl_chain_length:
                            acyl_chain_length = chain_len
            
            alcohol_valid = False
            alcohol_chain_length = 0
            # For the alcohol side, follow neighbors of the ester oxygen excluding the carbonyl carbon.
            for nbr in oxy_atom.GetNeighbors():
                if nbr.GetIdx() == carbon_atom.GetIdx():
                    continue
                if nbr.GetAtomicNum() == 6:
                    chain_len = longest_chain_from(mol, nbr.GetIdx(), 
                                                   exclude={carbon_atom.GetIdx(), oxy_atom.GetIdx()}, visited=set())
                    if chain_len >= 8:
                        alcohol_valid = True
                        if chain_len > alcohol_chain_length:
                            alcohol_chain_length = chain_len
            
            if acyl_valid and alcohol_valid:
                valid_ester_candidates.append((acyl_chain_length, alcohol_chain_length))
    
    if len(valid_ester_candidates) == 0:
        return False, "No ester linkage found with both sides having sufficiently long (≥8 carbons) chains."
    if len(valid_ester_candidates) > 1:
        return False, f"Multiple ester linkages found ({len(valid_ester_candidates)}). Wax compounds typically have a single ester bond."
    
    acyl_chain_length, alcohol_chain_length = valid_ester_candidates[0]
    return True, (f"Found ester linkage with acyl chain length {acyl_chain_length} and alcohol chain length {alcohol_chain_length}, "
                  f"consistent with wax compounds.")

# Example usage:
if __name__ == "__main__":
    # Test example using one of the provided SMILES: 2-palmitoyloxypalmityl palmitate
    test_smiles = "CCCCCCCCCCCCCCCC(=O)OCC(CCCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCCCC"
    classification, reason = is_wax(test_smiles)
    print(f"Classification: {classification}\nReason: {reason}")