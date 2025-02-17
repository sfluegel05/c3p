"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
"""
Classifies: trans-2-enoyl-CoA
An unsaturated fatty acyl-CoA that results from the formal condensation of the thiol group 
of coenzyme A with the carboxy group of any 2,3-trans-enoic acid.

This version first forces stereochemistry, then looks for a CoA core fragment while ensuring that
the CoA fragment has the expected overall formal charge (e.g. as in the neutralized form). It then 
identifies a thioester linkage in which the acyl side exhibits a double bond between the alpha 
and beta carbons that is marked as trans (E). Finally, it inspects the acyl chain length.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_trans_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a trans-2-enoyl-CoA based on its SMILES string.
    
    Steps:
      1. Validate the SMILES and assign stereochemistry.
      2. Look for a characteristic CoA core fragment. We use a SMARTS pattern that
         matches the key portion of the CoA backbone (ignoring minor charge differences)
         and then verify that the overall formal charge on the matched atoms is zero.
      3. Look for a thioester linkage defined as [CX3](=O)[SX2].
      4. For each thioester match, verify:
           • The carbonyl carbon is attached to exactly one carbon (the “alpha” carbon)
           • That alpha carbon shows a double bond to a beta carbon.
           • The double bond stereochemistry is trans (E).
           • The acyl chain is minimally extended (to rule out overly short chains).
         Also, check that the thioester sulfur is part of the CoA region.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is identified as trans-2-enoyl-CoA, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Force assignment of stereochemistry
    Chem.AssignStereochemistry(mol, force=True)
    
    # Look for a CoA core fragment.
    # We base the pattern on a conserved portion of the CoA backbone.
    # This pattern is intentionally less strict regarding phosphate charges.
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP")
    coa_matches = mol.GetSubstructMatches(coa_pattern)
    valid_coa_found = False
    for match in coa_matches:
        # Sum up the formal charges in this fragment
        if sum(mol.GetAtomWithIdx(idx).GetFormalCharge() for idx in match) == 0:
            valid_coa_found = True
            coa_match_indices = set(match)
            break
    if not valid_coa_found:
        return False, "Missing characteristic CoA moiety or improper protonation state"

    # Look for a thioester linkage defined as [CX3](=O)[SX2]
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Missing thioester (acyl-CoA linkage) moiety"
    
    # For each thioester match, examine the acyl chain
    for match in thioester_matches:
        # match tuple: index 0 is the carbonyl carbon, index 1 is the sulfur atom.
        carbonyl_idx = match[0]
        sulfur_idx = match[1]
        
        # Ensure that the thioester sulfur is part of the CoA core we found.
        if sulfur_idx not in coa_match_indices:
            # The thioester must connect the acyl chain to the CoA moiety.
            continue
        
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # Identify the "alpha" carbon attached to the carbonyl (exclude the sulfur).
        alpha_candidates = [nbr.GetIdx() for nbr in carbonyl_atom.GetNeighbors()
                            if nbr.GetIdx() != sulfur_idx and nbr.GetSymbol() == "C"]
        if not alpha_candidates:
            continue  # No valid alpha carbon found in this match.
        alpha_idx = alpha_candidates[0]
        alpha_atom = mol.GetAtomWithIdx(alpha_idx)
        
        # Look for a double bond emanating from the alpha carbon.
        acyl_double_found = False
        for bond in alpha_atom.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                beta_atom = bond.GetOtherAtom(alpha_atom)
                if beta_atom.GetSymbol() != "C":
                    continue
                # Check that the double bond is assigned as trans (E)
                stereo = bond.GetStereo()
                if stereo != Chem.rdchem.BondStereo.STEREOE:
                    return False, "Double bond adjacent to thioester is not in trans (E) configuration"
                
                # Extra check: Ensure that the alpha carbon is not overbranched.
                alpha_c_neighbors = [nbr for nbr in alpha_atom.GetNeighbors() if nbr.GetSymbol() == "C"]
                if len(alpha_c_neighbors) < 2:
                    return False, "Acyl chain too short: insufficient carbons attached to alpha carbon"
                
                # Walk along the acyl chain from the beta carbon and count carbon atoms.
                chain_length = 1  # count the beta carbon
                current_atom = beta_atom
                previous_idx = alpha_idx
                
                # Traverse up to 5 bonds to ensure a minimally long chain.
                for _ in range(5):
                    next_neighbors = [nbr for nbr in current_atom.GetNeighbors()
                                      if nbr.GetIdx() != previous_idx and nbr.GetSymbol() == "C"]
                    if not next_neighbors:
                        break
                    # follow the first neighbor (linear chain assumption)
                    previous_idx = current_atom.GetIdx()
                    current_atom = next_neighbors[0]
                    chain_length += 1
                
                if chain_length < 2:
                    return False, "Acyl chain appears too short"
                    
                acyl_double_found = True
                return True, ("Found thioester with an acyl chain exhibiting a trans (E) double bond at "
                              "the 2-enoyl position and a correctly protonated CoA moiety")
        # End bond loop on the alpha_atom
    # End thioester matches loop
    
    return False, "No appropriate trans-2-enoyl double bond detected in an acyl chain connected to CoA"

# For testing, you can run:
# if __name__ == "__main__":
#     test_smiles = "CCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
#     valid, reason = is_trans_2_enoyl_CoA(test_smiles)
#     print(valid, reason)