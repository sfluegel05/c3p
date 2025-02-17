"""
Classifies: CHEBI:73702 wax
"""
#!/usr/bin/env python3
"""
Classifies: Wax
Definition: A chemical substance that is an organic compound or mixture of compounds 
           that is composed of long-chain molecules and is malleable at ambient temperatures.
Approach:
  - Preliminary checks: must be organic (contain carbon),
    have a molecular weight ≥300 Da and at least 5 rotatable bonds.
  - Instead of using a SMARTS to simply locate an ester group,
    we iterate over all bonds to detect an ester linkage. An ester linkage is identified as:
      R–C(=O)–O–R'
    In our implementation, we look for a bond between an oxygen and a carbon where:
      • The carbon (ester carbon) is double-bonded to some other oxygen (the carbonyl oxygen).
      • Then we “cut” the molecule at the ester bond to consider separately the acyl (fatty acid) side
        and the alcoholic (fatty alcohol) side.
  - For each side, we perform a DFS over contiguous carbons (atomic num 6) to determine the maximum length
    of a linear carbon chain.
  - If we find at least one ester linkage where both sides yield a chain length of at least 8 carbons,
    we classify the molecule as a wax.
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
      exclude (set): Atom indices that we are not allowed to visit.
      visited (set): Atom indices already visited in the current path.
      
    Returns:
      int: The maximum number of carbons found in a contiguous chain (including the starting carbon).
    """
    visited.add(current_idx)
    max_length = 1  # count current carbon
    atom = mol.GetAtomWithIdx(current_idx)
    for nbr in atom.GetNeighbors():
        nbr_idx = nbr.GetIdx()
        # Do not cross into excluded atoms; avoid cycles.
        if nbr_idx in visited or nbr_idx in exclude:
            continue
        # Only consider carbons
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
      - Must be organic (contain carbon atoms).
      - Must have molecular weight ≥300 Da and at least 5 rotatable bonds.
      - Must contain at least one ester linkage (R–C(=O)–O–R') where, after "cutting" the molecule at the ester bond,
        both the acyl (fatty acid) part and the alcohol (fatty alcohol) part yield a longest contiguous carbon 
        chain of at least 8 atoms.
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if molecule qualifies as a wax, False otherwise.
      str: Detailed reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization failed: {e}"
    
    # Check for organic molecule (must have at least one carbon)
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Not organic (no carbon atoms found)."
    
    # Check molecular weight and rotatable bonds
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight ({mol_wt:.1f} Da) too low for a typical wax."
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, f"Too few rotatable bonds ({n_rotatable}) for long-chain wax characteristics."
    
    # Now search for ester linkages: R-C(=O)-O-R'
    # Iterate over all bonds and look for a bond between an oxygen and a carbon.
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Identify potential ester bond: one atom is oxygen, the other is carbon.
        if (a1.GetAtomicNum() == 8 and a2.GetAtomicNum() == 6) or (a2.GetAtomicNum() == 8 and a1.GetAtomicNum() == 6):
            # Assign roles: carbon should be the carbonyl carbon if it is double-bonded to another oxygen.
            if a1.GetAtomicNum() == 6:
                carbon_atom = a1
                oxy_atom = a2
            else:
                carbon_atom = a2
                oxy_atom = a1
            
            # Check that carbon_atom is a carbonyl carbon: it must be double-bonded to some oxygen 
            # (other than the current oxy_atom) to qualify as an ester carbon.
            carbonyl_found = False
            for nbr in carbon_atom.GetNeighbors():
                if nbr.GetIdx() == oxy_atom.GetIdx():
                    continue
                # Look for an oxygen with a double bond.
                bond_order = mol.GetBondBetweenAtoms(carbon_atom.GetIdx(), nbr.GetIdx()).GetBondTypeAsDouble()
                if nbr.GetAtomicNum() == 8 and bond_order >= 2.0:
                    carbonyl_found = True
                    break
            if not carbonyl_found:
                continue  # not an ester linkage
            
            # Now we have recognized an ester linkage.
            # For the acyl side (fatty acid), we follow the carbonyl carbon’s neighbors 
            # EXCEPT the ester oxygen (and the carbonyl oxygen that is double-bonded).
            acyl_valid = False
            acyl_chain_length = 0
            for nbr in carbon_atom.GetNeighbors():
                if nbr.GetIdx() == oxy_atom.GetIdx():
                    continue
                # Exclude the carbonyl oxygen (which is double-bonded, if any)
                bond_tmp = mol.GetBondBetweenAtoms(carbon_atom.GetIdx(), nbr.GetIdx())
                if nbr.GetAtomicNum() == 8 and bond_tmp.GetBondTypeAsDouble() >= 2.0:
                    continue
                if nbr.GetAtomicNum() == 6:
                    chain_len = longest_chain_from(mol, nbr.GetIdx(), exclude={oxy_atom.GetIdx(), carbon_atom.GetIdx()}, visited=set())
                    if chain_len >= 8:
                        acyl_valid = True
                        if chain_len > acyl_chain_length:
                            acyl_chain_length = chain_len
            # For the alcoholic side (fatty alcohol), we follow the ester oxygen's neighbors excluding the carbon.
            alcohol_valid = False
            alcohol_chain_length = 0
            for nbr in oxy_atom.GetNeighbors():
                if nbr.GetIdx() == carbon_atom.GetIdx():
                    continue
                if nbr.GetAtomicNum() == 6:
                    chain_len = longest_chain_from(mol, nbr.GetIdx(), exclude={carbon_atom.GetIdx(), oxy_atom.GetIdx()}, visited=set())
                    if chain_len >= 8:
                        alcohol_valid = True
                        if chain_len > alcohol_chain_length:
                            alcohol_chain_length = chain_len
            if acyl_valid and alcohol_valid:
                return True, (f"Found ester linkage with acyl chain length {acyl_chain_length} "
                              f"and alcohol chain length {alcohol_chain_length}, consistent with wax compounds.")
    
    return False, "No ester linkage found with both sides having sufficiently long (≥8 carbons) chains."

# Example usage:
if __name__ == "__main__":
    # Test with one of the provided SMILES: 2-palmitoyloxypalmityl palmitate
    test_smiles = "CCCCCCCCCCCCCCCC(=O)OCC(CCCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCCCC"
    classification, reason = is_wax(test_smiles)
    print(f"Classification: {classification}\nReason: {reason}")