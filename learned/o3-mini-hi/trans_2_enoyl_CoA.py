"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
#!/usr/bin/env python3
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
    
    The function performs the following checks:
      1. Validity of the SMILES.
      2. Presence of a CoA moiety (via a common CoA substructure fragment).
      3. Presence of a thioester linkage ([CX3](=O)[SX2]), which defines the acyl-CoA connection.
      4. The acyl chain (attached to the carbonyl carbon of the thioester) must include a
         double bond between the alpha (C2) and beta (C3) carbons with trans (E) configuration.
         
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is identified as trans-2-enoyl-CoA, False otherwise.
        str: Reason for classification.
    """
    # Parse the molecule from the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the Coenzyme A (CoA) moiety.
    # Here we use a characteristic CoA fragment: "SCCNC(=O)CCNC(=O)"
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing characteristic CoA moiety"
    
    # Look for a thioester moiety defined as [CX3](=O)[SX2]
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Missing thioester (acyl-CoA linkage) moiety"
    
    # Loop through each thioester match. The match returns a tuple of atom indices: (carbonyl, sulfur)
    for match in thioester_matches:
        carbonyl_idx = match[0]
        sulfur_idx = match[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Find the 'alpha' carbon: a carbon attached to the carbonyl that is not the sulfur.
        alpha_candidates = [nbr.GetIdx() for nbr in carbonyl_atom.GetNeighbors() 
                            if nbr.GetIdx() != sulfur_idx and nbr.GetSymbol() == "C"]
        if not alpha_candidates:
            continue  # Skip if no candidate found
        
        alpha_idx = alpha_candidates[0]
        alpha_atom = mol.GetAtomWithIdx(alpha_idx)
        
        # Look for a double bond from the alpha carbon to a beta carbon.
        # In a 2-enoyl chain the double bond should be directly connected to the alpha carbon.
        for bond in alpha_atom.GetBonds():
            # We want a double bond.
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                # Identify the beta carbon (the other atom in the double bond, not the alpha carbon)
                beta_atom = bond.GetOtherAtom(alpha_atom)
                # Sanity check: beta_atom should be a carbon.
                if beta_atom.GetSymbol() != "C":
                    continue
                # Check the stereochemistry of the double bond.
                # RDKit assigns stereochemistry (STEREOE for trans, STEREOZ for cis).
                stereo = bond.GetStereo()
                if stereo == Chem.rdchem.BondStereo.STEREOE:
                    return True, "Found thioester with an acyl chain exhibiting a trans (E) double bond at the 2-enoyl position and a CoA moiety"
                else:
                    # If a double bond is found but not marked as trans, that is a potential reason.
                    return False, "Double bond found adjacent to the thioester is not in trans (E) configuration"
    
    return False, "No appropriate trans-2-enoyl double bond detected in the acyl chain"

# For testing purposes (uncomment the lines below):
# if __name__ == "__main__":
#     test_smiles = "CCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC..."  # example from list
#     valid, reason = is_trans_2_enoyl_CoA(test_smiles)
#     print(valid, reason)