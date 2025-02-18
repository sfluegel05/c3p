"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
#!/usr/bin/env python3
"""
Classifies: 3-oxo-fatty acyl-CoA
An oxo fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A 
with the carboxy group of any 3-oxo-fatty acid.
The idea is to find a fatty acyl thioester fragment of the form:
   R–C(=O)–C(R’)–C(=O)S– (where the C(R’) is the 3-oxo carbon)
and then to ensure that the C(R’) is not additionally substituted by an extra carboxyl group.
We also require the presence of the typical CoA adenine fragment.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.
    
    Criteria:
      1. The molecule must be valid and not contain explicit deprotonated oxygens (i.e. [O-])
         so that we only classify neutral CoA derivatives.
      2. The molecule should contain a 3-oxo fatty acyl thioester fragment. We define a SMARTS pattern
         for a fragment of the form R–C(=O)–C(R’)–C(=O)S, where the central carbon (C(R’)) must be
         aliphatic. (The previous version’s pattern was too loose.)
      3. In any such match the middle carbon must not have an extra carboxyl substituent.
      4. The molecule must also contain a Coenzyme A signature (here we require the adenine fragment).
    
    Args:
      smiles (str): SMILES string representing the molecule.
    
    Returns:
      bool: True if the molecule is classified as a 3-oxo-fatty acyl-CoA, False otherwise.
      str: Explanation for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Exclude molecules that have explicit deprotonated oxygens (e.g. [O-])
    if "[O-]" in smiles:
        return False, "Contains deprotonated oxygens; expected neutral CoA"
    
    # First, check that the molecule contains a CoA adenine fragment.
    coa_pattern = Chem.MolFromSmarts("n1cnc2ncnc(c12)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Does not contain a Coenzyme A adenine fragment"
    
    # Define a refined SMARTS for the 3-oxo fatty acyl thioester fragment.
    # This pattern looks for a fragment of the form: [#6]-C(=O)-[#6X4]-C(=O)S
    # where the atom tagged as [#6X4] (the middle carbon) should be aliphatic.
    frag_smarts = "[#6:1]-C(=O)-[#6X4:2]-C(=O)S"
    frag_pattern = Chem.MolFromSmarts(frag_smarts)
    matches = mol.GetSubstructMatches(frag_pattern)
    if not matches:
        return False, "Does not contain a 3-oxo fatty acyl thioester substructure"
    
    # For each match, we want to ensure that the middle carbon (tag 2) is not substituted
    # with an extra carboxyl group (e.g. -C(=O)O). In a proper 3-oxo fatty acyl moiety the middle
    # carbon is only part of the two carbonyls.
    #
    # We'll get the middle atom from the match and look at its neighbors (excluding those that
    # are part of the fragment match). If one of them is a carbon that is in a carbonyl (i.e. has a 
    # double bond from that carbon to oxygen), we consider that an extra carboxyl.
    def middle_has_extra_carboxyl(mid_atom, exclude_idxs):
        # Iterate over neighbors of the middle atom that are not in the exclude list.
        for nbr in mid_atom.GetNeighbors():
            if nbr.GetIdx() in exclude_idxs:
                continue
            # Check: if the neighbor is a carbon that itself is doubly bonded to oxygen,
            # then count that as a carboxyl substituent.
            if nbr.GetAtomicNum() == 6:
                for bond in nbr.GetBonds():
                    # bond between nbr and an oxygen that is a double bond
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        other = bond.GetOtherAtom(nbr)
                        if other.GetAtomicNum() == 8:
                            return True
        return False

    valid_match_found = False
    for match in matches:
        # match is a tuple of atom indices corresponding to the SMARTS pattern:
        # index0: the first carbon (R group), index1: its carbonyl carbon, 
        # index2: the middle (aliphatic) carbon, index3: the carbonyl carbon of the thioester.
        mid_idx = match[2]
        # Exclude the two atoms that are part of the acyl fragment (index1 and index3)
        exclude = {match[1], match[3]}
        mid_atom = mol.GetAtomWithIdx(mid_idx)
        if middle_has_extra_carboxyl(mid_atom, exclude):
            continue
        else:
            valid_match_found = True
            break

    if not valid_match_found:
        return False, "Contains an acyl thioester-like fragment but with extra carboxyl substitution (e.g. malonyl-type) rather than a 3-oxo fatty acyl group"
    
    # Passed all checks:
    return True, "Contains a valid 3-oxo fatty acyl thioester linked to a Coenzyme A moiety"

# Example usage:
if __name__ == "__main__":
    # Example SMILES for 3-oxoadipyl-CoA:
    smiles_example = ("CC(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1"
                      "OP(O)(O)=O)n1cnc2c(N)ncnc12)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)"
                      "CC(=O)CCC(O)=O")
    result, reason = is_3_oxo_fatty_acyl_CoA(smiles_example)
    print(result, reason)