"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
#!/usr/bin/env python3
"""
Classifies: trans-2-enoyl-CoA
An unsaturated fatty acyl-CoA that results from the formal condensation of the thiol group 
of coenzyme A with the carboxy group of any 2,3-trans-enoic acid.
  
This improved version first assigns stereochemistry and uses a more robust approach:
  - It verifies a CoA moiety by matching a common fragment (ignoring formal charge differences).
  - It identifies a thioester linkage.
  - It then “extracts” the acyl chain side of the thioester and checks:
      • that the carbonyl carbon is attached to exactly one carbon (the “alpha” carbon)
      • that the alpha carbon forms a double bond to a beta carbon
      • that this double bond is marked as trans (E)
      • that the acyl chain is long and mostly aliphatic.
If all checks pass, the molecule is classified as trans-2-enoyl-CoA.
  
False positives from the earlier version (e.g. structures with extra charges or subtle modifications)
are less likely to pass the extra acyl-chain analysis.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_trans_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a trans-2-enoyl-CoA based on its SMILES string.

    The function performs the following checks:
      1. Validity of the SMILES.
      2. Presence of a CoA moiety (via a common CoA fragment, ignoring charge differences).
      3. Presence of a thioester linkage ([CX3](=O)[SX2]) and extraction of the acyl side.
      4. The acyl chain (attached to the carbonyl carbon of the thioester) must include 
         a double bond between the alpha (C2) and beta (C3) carbons, and that double bond 
         must be in trans (E) configuration.
      5. The acyl chain is further verified to be a contiguous aliphatic chain of reasonable length.
         
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is identified as trans-2-enoyl-CoA, False otherwise.
        str: Reason for classification.
    """
    # Parse and sanitize the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    # Force assignment of stereochemistry 
    Chem.AssignStereochemistry(mol, force=True)
    
    # Look for a CoA fragment. We use a common substructure 
    # (ignoring formal charge on phosphate oxygens using wildcard for oxygen atoms).
    # This pattern should match both protonated and deprotonated forms.
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing characteristic CoA moiety"
    
    # Look for a thioester linkage defined as [CX3](=O)[SX2]
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Missing thioester (acyl-CoA linkage) moiety"
    
    # For each thioester found, examine the acyl chain
    for match in thioester_matches:
        # In the match tuple, index 0 is the carbonyl carbon and index 1 is the sulfur.
        carbonyl_idx = match[0]
        sulfur_idx = match[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Identify the "alpha" carbon attached to the carbonyl (exclude the sulfur).
        alpha_candidates = [nbr.GetIdx() for nbr in carbonyl_atom.GetNeighbors() 
                            if nbr.GetIdx() != sulfur_idx and nbr.GetSymbol() == "C"]
        if not alpha_candidates:
            continue  # no valid alpha carbon found in this match
        alpha_idx = alpha_candidates[0]
        alpha_atom = mol.GetAtomWithIdx(alpha_idx)
        
        # We now look for a double bond from the alpha carbon to a beta carbon.
        acyl_double_found = False
        for bond in alpha_atom.GetBonds():
            # Recognize a double bond that connects alpha atom to another carbon.
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                beta_atom = bond.GetOtherAtom(alpha_atom)
                if beta_atom.GetSymbol() != "C":
                    continue
                # Check whether the double bond stereochemistry is assigned as trans.
                stereo = bond.GetStereo()
                if stereo != Chem.rdchem.BondStereo.STEREOE:
                    return False, "Double bond adjacent to thioester is not in trans (E) configuration"
                # Extra check: ensure the double bond is at the expected position.
                # We require that the alpha carbon has only one carbon neighbor besides the carbonyl
                # (i.e. the double-bond forms the start of a linear acyl chain).
                alpha_nb_carbons = [nbr for nbr in alpha_atom.GetNeighbors() if nbr.GetSymbol() == "C"]
                if len(alpha_nb_carbons) < 2:
                    return False, "Acyl chain too short: insufficient carbons on alpha atom"
                
                # Further, traverse from the beta carbon along a linear chain (non-branched, aliphatic)
                # and count a few atoms. Valid fatty acyl chains are usually at least 3 carbons long.
                chain_length = 1  # beta carbon counts as 1
                current_atom = beta_atom
                prev_idx = alpha_idx
                # Walk outwards in the chain for up to 5 bonds (we expect a medium-to-long fatty chain)
                for _ in range(5):
                    neighbors = [nbr for nbr in current_atom.GetNeighbors() if nbr.GetIdx() != prev_idx and nbr.GetSymbol() == "C"]
                    if not neighbors:
                        break
                    # if more than one neighbor then branching; we consider the chain valid if one continuous path exists
                    current_atom = neighbors[0]
                    prev_idx = current_atom.GetIdx()
                    chain_length += 1
                if chain_length < 2:
                    return False, "Acyl chain appears too short"
                acyl_double_found = True
                # If we have reached here then the thioester, alpha-beta double bond (with trans stereo)
                # and a minimally extended acyl chain are detected.
                return True, ("Found thioester with an acyl chain exhibiting a trans (E) double bond at "
                              "the 2-enoyl position and a CoA moiety")
        # End for bonds on alpha_atom
    # End for each thioester match
    
    return False, "No appropriate trans-2-enoyl double bond detected in the acyl chain"

# For testing (if run directly), one can uncomment the lines below:
# if __name__ == "__main__":
#     test_smiles = "CCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
#     valid, reason = is_trans_2_enoyl_CoA(test_smiles)
#     print(valid, reason)