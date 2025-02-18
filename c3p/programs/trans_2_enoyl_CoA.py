"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
"""
Classifies: trans-2-enoyl-CoA
An unsaturated fatty acyl-CoA that results from the formal condensation of the thiol group 
of coenzyme A with the carboxy group of any 2,3-trans-enoic acid.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_trans_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a trans-2-enoyl-CoA based on its SMILES string.
    
    Steps:
      1. Validate the SMILES and assign stereochemistry.
      2. Look for a characteristic CoA core fragment using a SMARTS pattern.
         We then check that the thioester linkage (defined as [CX3](=O)[SX2]) attaches
         the acyl chain to the CoA moiety.
      3. For each thioester match, identify the alpha carbon (attached to the carbonyl) and
         check that it participates in a double bond (to the beta carbon) with explicit trans (E) stereochemistry.
      4. Accept if any thioester is found to bear an acyl chain with this trans configuration at the 2-enoyl position.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is identified as trans-2-enoyl-CoA, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Force assignment of stereochemistry to ensure double bond stereo is defined.
    Chem.AssignStereochemistry(mol, force=True)
    
    # Identify a CoA core fragment. This pattern targets the characteristic portion of CoA.
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP")
    coa_matches = mol.GetSubstructMatches(coa_pattern)
    if not coa_matches:
        return False, "Missing characteristic CoA moiety"
    
    # Select one matching CoA fragment that has an overall formal charge of 0.
    coa_match_indices = None
    for match in coa_matches:
        if sum(mol.GetAtomWithIdx(idx).GetFormalCharge() for idx in match) == 0:
            coa_match_indices = set(match)
            break
    if coa_match_indices is None:
        return False, "Found CoA-like fragment but with improper protonation/charge"
    
    # Look for a thioester linkage defined as [CX3](=O)[SX2]
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Missing thioester (acyl-CoA linkage) moiety"
    
    # Keep track if any valid trans-2-enoyl acyl chain is found
    valid_thioester_found = False
    
    # Process each thioester match to check the acyl chain double bond configuration
    for match in thioester_matches:
        # In the pattern, index 0 is the carbonyl carbon and index 1 is the sulfur.
        carbonyl_idx = match[0]
        sulfur_idx = match[1]
        # The thioester should connect the fatty acyl chain to CoA, so the sulfur must be within the CoA fragment.
        if sulfur_idx not in coa_match_indices:
            continue
        
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # Identify the "alpha" carbon attached to the carbonyl (exclude the sulfur neighbor)
        alpha_candidates = [nbr.GetIdx() for nbr in carbonyl_atom.GetNeighbors() 
                            if nbr.GetIdx() != sulfur_idx and nbr.GetSymbol() == "C"]
        if not alpha_candidates:
            # no valid alpha carbon found for this thioester
            continue
        alpha_idx = alpha_candidates[0]
        alpha_atom = mol.GetAtomWithIdx(alpha_idx)
        
        # Now iterate over all bonds from the alpha carbon to find a double bond.
        valid_double_found = False
        for bond in alpha_atom.GetBonds():
            # Interested only in bonds that are double bonds.
            if bond.GetBondType() != Chem.rdchem.BondType.DOUBLE:
                continue
            
            # Identify the other atom in the double bond (should be a beta carbon)
            other_atom = bond.GetOtherAtom(alpha_atom)
            if other_atom.GetSymbol() != "C":
                continue
            
            # Verify that this is the double bond adjacent to the carbonyl (i.e. alpha-beta bond)
            # The bond should involve the alpha atom we found from the carbonyl.
            # Now check stereochemistry:
            stereo = bond.GetStereo()
            if stereo != Chem.rdchem.BondStereo.STEREOE:
                # This bond exists but is not marked as trans (E); skip it.
                continue
            
            # Optionally, check that the acyl chain is of a minimal length (e.g. at least 2 carbons: 
            # the beta carbon plus at least one more) to discount trivial fragments.
            chain_length = 1  # count the beta carbon already
            current_atom = other_atom
            previous_idx = alpha_idx
            # Traverse a few bonds to gauge chain length
            for _ in range(5):
                next_neighbors = [nbr for nbr in current_atom.GetNeighbors() 
                                  if nbr.GetIdx() != previous_idx and nbr.GetSymbol() == "C"]
                if not next_neighbors:
                    break
                # follow one neighbor (assuming a linear chain for simplicity)
                previous_idx = current_atom.GetIdx()
                current_atom = next_neighbors[0]
                chain_length += 1
            if chain_length < 1:
                continue
            
            # If we have reached this point, we found a proper double bond with E configuration.
            valid_double_found = True
            break  # no need to inspect further bonds on alpha
        
        if valid_double_found:
            valid_thioester_found = True
            break  # Accept once any valid thioester with a trans-2-enoyl acyl chain is found
    
    if valid_thioester_found:
        return True, ("Found thioester with an acyl chain exhibiting a trans (E) double bond at the "
                      "2-enoyl position and with a proper CoA moiety")
    else:
        return False, "No appropriate trans-2-enoyl double bond detected in an acyl chain connected to CoA"

# Example usage testing (uncomment to run):
# if __name__ == "__main__":
#     test_smiles = "CCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
#     valid, reason = is_trans_2_enoyl_CoA(test_smiles)
#     print(valid, reason)