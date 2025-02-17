"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
"""
Classifies: trans-2-enoyl-CoA
An unsaturated fatty acyl-CoA that results from the formal condensation of the thiol group 
of coenzyme A with the carboxy group of any 2,3-trans-enoic acid.
 
The revised approach:
  1. Validate the SMILES and force assignment of stereochemistry.
  2. Identify a CoA motif using a looser SMARTS pattern.
  3. Find thioester linkages (the [CX3](=O)[SX2] motif) that connect the acyl chain to CoA.
  4. For each thioester, from the carbonyl carbon identify the acyl chain’s alpha carbon (i.e. the one not attached to sulfur).
  5. From the alpha carbon, attempt to find a double bond (alpha-beta bond) that is directly attached.
     We then verify that the bond is defined with trans (E) stereochemistry.
  6. Also (optionally) traverse a short distance along the acyl chain to ensure that the chain is of a minimal length.
 
If a thioester with an acyl chain that has a proper trans double bond in the 2-enoyl position is found, we return True.
Otherwise we return False with an explanation.
 
Note: This is a heuristic procedure and may require further tuning.
"""
 
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
 
def is_trans_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a trans-2-enoyl-CoA based on its SMILES string.
 
    Args:
        smiles (str): SMILES string of the molecule.
 
    Returns:
        bool: True if the molecule is identified as a trans-2-enoyl-CoA, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
 
    # Force assignment of stereochemistry (this will help for bonds defined with / and \)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
 
    # Identify a CoA fragment using a less-restrictive pattern.
    # Here we look for a substructure commonly present in CoA (many examples have 'SCCNC(=O)CCNC(=O)')
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing characteristic CoA motif"
 
    # Find the thioester linkage: a carbonyl attached to a sulfur.
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Missing thioester (acyl-CoA linkage) moiety"
 
    # Loop over each thioester match to check for a valid trans-2-enoyl chain.
    # In the match, index0 is carbonyl carbon and index1 is sulfur.
    for match in thioester_matches:
        carbonyl_idx = match[0]
        sulfur_idx = match[1]
 
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # Find alpha carbon: a neighbor of the carbonyl that is not the sulfur and is a carbon.
        alpha_candidates = [nbr.GetIdx() for nbr in carbonyl_atom.GetNeighbors() 
                            if nbr.GetIdx() != sulfur_idx and nbr.GetAtomicNum() == 6]
        if not alpha_candidates:
            continue
        alpha_idx = alpha_candidates[0]
        alpha_atom = mol.GetAtomWithIdx(alpha_idx)
 
        # Now, look for a double bond leaving the alpha carbon.
        valid_double_found = False
        for bond in alpha_atom.GetBonds():
            # We are interested only in bonds that are double bonds.
            if bond.GetBondType() != Chem.rdchem.BondType.DOUBLE:
                continue
 
            # The double bond should be between the alpha carbon and the beta carbon.
            beta_atom = bond.GetOtherAtom(alpha_atom)
            if beta_atom.GetAtomicNum() != 6:
                continue
 
            # For a 2-enoyl chain, the alpha carbon should have a single (unique) double bond to a beta carbon.
            # Check the stereochemistry of this double bond.
            # Note: The stereo information should have been assigned above.
            bond_stereo = bond.GetStereo()
            if bond_stereo != Chem.rdchem.BondStereo.STEREOE:
                # Not a trans bond so skip.
                continue
 
            # Optionally, check that the acyl chain is at least a minimal length.
            chain_length = 1  # count beta carbon
            current_atom = beta_atom
            previous_idx = alpha_idx
            # Traverse up to 4 additional bonds to ensure chain length.
            for _ in range(4):
                next_neighbors = [nbr for nbr in current_atom.GetNeighbors() 
                                  if nbr.GetIdx() != previous_idx and nbr.GetAtomicNum() == 6]
                if not next_neighbors:
                    break
                previous_idx = current_atom.GetIdx()
                current_atom = next_neighbors[0]
                chain_length += 1
            if chain_length < 1:  # this check is trivial; in practice a real acyl chain will have several carbons
                continue
 
            # We have found a double bond from the alpha carbon with E stereochemistry.
            valid_double_found = True
            break  # no need to check further bonds from this alpha atom
 
        if valid_double_found:
            return True, ("Found thioester linkage with an acyl chain exhibiting a trans (E) double bond "
                          "at the 2-enoyl position and proper CoA signature")
 
    return False, "No appropriate trans-2-enoyl double bond detected in an acyl chain connected to CoA"
 
# Example usage (uncomment to run tests):
# if __name__ == "__main__":
#     test_smiles = "CCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
#     valid, reason = is_trans_2_enoyl_CoA(test_smiles)
#     print(valid, reason)