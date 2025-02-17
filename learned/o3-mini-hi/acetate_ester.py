"""
Classifies: CHEBI:47622 acetate ester
"""
"""
Classifies: Acetate ester – Any carboxylic ester where the acyl (carboxylic acid) component is acetic acid.
That is, the molecule must contain an ester group whose acyl part is CH3C(=O)–.
We identify this by looping over all carbonyl carbons and checking if:
  • There is one double-bonded oxygen (i.e. the C=O oxygen).
  • There is one single-bonded oxygen (i.e. the ester oxygen) that is not –OH.
  • The remaining substituent on the carbonyl carbon is a methyl group (CH3).
"""

from rdkit import Chem

def is_acetate_ester(smiles: str):
    """
    Determines if a molecule is an acetate ester based on its SMILES string.
    An acetate ester is any carboxylic ester with acyl component CH3C(=O)-.
    This implementation loops over all carbons that are part of a carbonyl group,
    then checks if one of their heavy-atom neighbors is an oxygen (ester oxygen)
    and the remaining neighbor is a CH3 group.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule contains at least one validated acetate ester moiety, False otherwise.
        str: Explanation for the decision.
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    # Ensure hydrogen counts are updated
    mol.UpdatePropertyCache(strict=False)
    
    acetate_found = False
    
    # Iterate over all atoms looking for a carbon candidate for the acyl carbon.
    for atom in mol.GetAtoms():
        # Look at heavy atoms that are carbons only.
        if atom.GetAtomicNum() != 6:
            continue
        
        # In an ester acyl group, the carbonyl carbon (C=O) is sp2 and connected to three heavy atoms:
        # one oxygen by a double bond (C=O), one oxygen by a single bond (the ester oxygen),
        # and one carbon as the alkyl (methyl) group.
        # So we first require that the number of neighbors (heavy atoms) is exactly 3.
        if atom.GetDegree() != 3:
            continue

        neighbors = atom.GetNeighbors()
        
        dbl_oxygen = None
        ester_oxygen = None
        alkyl_group = None
        
        # Identify which neighbor is the double-bond oxygen.
        # We check bonds to see if one bond is double and the neighbor is oxygen.
        for nbr in neighbors:
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            # Check if the neighbor is oxygen
            if nbr.GetAtomicNum() == 8:
                if bond.GetBondTypeAsDouble() == 2.0:
                    dbl_oxygen = nbr
                elif bond.GetBondType() == Chem.BondType.SINGLE:
                    ester_oxygen = nbr
        
        # The remaining neighbor (if any) must be the alkyl part.
        # It is possible that the order of neighbors is not as expected so we try to identify it.
        if dbl_oxygen is not None and ester_oxygen is not None:
            for nbr in neighbors:
                # Skip the ones used as oxygens
                if nbr.GetIdx() in (dbl_oxygen.GetIdx(), ester_oxygen.GetIdx()):
                    continue
                # This candidate should be a carbon
                if nbr.GetAtomicNum() == 6:
                    alkyl_group = nbr
                    break
        
        # If any of the parts is missing, skip this atom.
        if dbl_oxygen is None or ester_oxygen is None or alkyl_group is None:
            continue
            
        # Check that the alkyl group is a methyl group:
        # In RDKit, a CH3 group should have a total hydrogen count of 3.
        if alkyl_group.GetTotalNumHs() != 3:
            continue

        # Ensure that the oxygen attached via a single bond (ester oxygen) is not a hydroxyl,
        # i.e. it should NOT have any attached hydrogen.
        # (In an ester, the oxygen should be bridging to another heavy atom.)
        if ester_oxygen.GetTotalNumHs() != 0:
            continue

        # If we reach here, we have found a candidate:
        # atom = carbonyl carbon with a C=O (dbl_oxygen), a single-bonded oxygen (ester_oxygen),
        # and an alkyl substituent which is CH3.
        acetate_found = True
        break

    if acetate_found:
        return True, "Contains validated acetate ester moiety (CH3C(=O)O-)"
    else:
        return False, "Does not contain a validated acetate ester moiety (CH3C(=O)O-)"

    
# Example usage:
if __name__ == "__main__":
    # Test some examples
    # methyl acetate should be classified as a valid acetate ester.
    test_smiles_list = [
        ("COC(=O)C", "methyl acetate"),
        ("CC(=O)OC", "methyl acetate alternative"),
        ("O=C(O)C", "acetic acid (not an ester)"),
    ]
    for smi, name in test_smiles_list:
        result, reason = is_acetate_ester(smi)
        print(name, smi, "->", result, reason)