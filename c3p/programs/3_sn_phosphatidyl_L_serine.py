"""
Classifies: CHEBI:11750 3-sn-phosphatidyl-L-serine
"""
"""
Classifies: 3-sn-phosphatidyl-L-serine 
Definition: 
  A 3-sn-glycerophosphoserine compound having acyl substituents at the 1- and 2-hydroxy positions.
  
Heuristic approach:
1. Parse the SMILES.
2. Check if the molecule contains a phosphoserine head group (using a SMARTS that matches “COP(=O)(O)OC[C](N)C(=O)O”).
3. Look for ester bonds where the acyl group (R-CO–) is attached to the glycerol backbone. 
   We identify ester bonds by matching the pattern "[CX3](=O)O". For each match, we retrieve the 
   carbonyl carbon and then follow its substituent (not the double-bonded oxygen) to judge if it 
   corresponds to an acyl chain that is “long” (here we require at least 6 carbon atoms along a chain).
4. The compound must have exactly two such acyl chains.
5. If any step fails, an appropriate message is returned.
  
Note: This is only one possible heuristic implementation.
"""

from rdkit import Chem

def is_3_sn_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a 3-sn-phosphatidyl-L-serine based on its SMILES string.
    The compound must have a glycerophosphoserine headgroup plus two acyl substituents at the
    1- and 2-hydroxy positions.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is 3-sn-phosphatidyl-L-serine, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for the presence of phosphorus (P) and nitrogen (N)
    # (since the phosphoserine headgroup must contain these atoms)
    hasP = any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms())
    hasN = any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms())
    if not hasP or not hasN:
        return False, "Missing key atoms (P and/or N) required for phosphoserine"

    # Define a SMARTS pattern for the phosphoserine head-group.
    # This pattern looks for an –O–P(=O)(O)O–C–C(N)C(=O)O fragment.
    # (Chirality markers are omitted for a general match.)
    ps_smarts = "CO[P](=O)(O)OC[C](N)C(=O)O"
    ps_pattern = Chem.MolFromSmarts(ps_smarts)
    if not ps_pattern:
        return False, "Internal error processing phosphoserine SMARTS"
    
    if not mol.HasSubstructMatch(ps_pattern):
        return False, "Phosphoserine head-group pattern not found"
        
    # Now, search for acyl ester bonds.
    # We look for ester bonds using the pattern "[CX3](=O)O".
    ester_smarts = "[CX3](=O)O"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester bonds found"
    
    # Define helper function to recursively compute longest chain length starting from a given atom.
    def longest_carbon_chain(mol, atom_idx, visited):
        atom = mol.GetAtomWithIdx(atom_idx)
        # Only count carbon atoms
        if atom.GetAtomicNum() != 6:
            return 0
        max_length = 1
        visited.add(atom_idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in visited:
                continue
            # Proceed only if neighbor is carbon 
            if nbr.GetAtomicNum() == 6:
                chain_length = 1 + longest_carbon_chain(mol, nbr_idx, visited.copy())
                if chain_length > max_length:
                    max_length = chain_length
        return max_length

    # For each ester match, determine if it likely represents an acyl substituent
    # attached to a glycerol backbone (the acyl group will be the fragment attached to the carbonyl carbon).
    acyl_count = 0
    for match in ester_matches:
        # The SMARTS "[CX3](=O)O" returns three atoms:
        # match[0] = carbonyl carbon, match[1] = double-bond oxygen, match[2] = ester oxygen.
        carbonyl_idx = match[0]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # Find the neighbor of the carbonyl carbon that is not the double-bond oxygen.
        acyl_neighbor = None
        for nbr in carbonyl_atom.GetNeighbors():
            if nbr.GetIdx() == match[1]:
                continue  # skip double-bonded oxygen
            # In an ester, the carbonyl carbon is attached to the acyl chain (R) and the ester oxygen.
            # We want the substituent R (the acyl chain).
            acyl_neighbor = nbr
            break
        if acyl_neighbor is None:
            continue
        # Compute the length of the carbon chain (number of carbons) in the acyl substituent.
        chain_length = longest_carbon_chain(mol, acyl_neighbor.GetIdx(), set())
        # We require a minimum chain length (here arbitrarily using 6 carbons) to qualify as an acyl chain.
        if chain_length >= 6:
            acyl_count += 1

    if acyl_count != 2:
        return False, f"Expected 2 acyl groups with sufficiently long chains; found {acyl_count}"
    
    # If all tests pass, it is classified as a 3-sn-phosphatidyl-L-serine.
    return True, "Molecule contains a phosphoserine head-group and two acyl substituents at the 1- and 2-hydroxy positions"

# Example usage (you can remove or comment out the lines below when using this module):
if __name__ == "__main__":
    test_smiles = "CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)OC[C@H](N)C(O)=O)OC(=O)CCCCCCCCCCCCCCCCC"
    result, reason = is_3_sn_phosphatidyl_L_serine(test_smiles)
    print(result, reason)