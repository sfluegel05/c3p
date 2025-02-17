"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
#!/usr/bin/env python3
"""
Classifies: 3-oxo-fatty acyl-CoA
An oxo fatty acyl-CoA results from the condensation of the thiol group of coenzyme A
with the carboxy group of any 3-oxo-fatty acid.
This program looks for a fatty acyl thioester fragment of the form:
   R–C(=O)–C(=O)S–
where the “middle” carbon (the one between the two carbonyls) should not be sterically overloaded.
Because the central carbon may be sp3 or sp2 (for example, when part of a ring)
we try two SMARTS patterns.
We also require the presence of a CoA adenine fragment.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.

    Criteria:
      1. Input must parse as a molecule and not contain deprotonated oxygens (i.e. [O-])
         so that only neutral CoA derivatives are accepted.
      2. The molecule should contain a 3-oxo fatty acyl thioester fragment. In a valid fragment,
         a pattern of the type: [#6]-C(=O)-C(=O)S is required. However, the “middle” carbon (the one
         between the two carbonyls) may be sp3 or sp2 (especially when it is embedded in alicyclic
         systems). We therefore try two patterns.
      3. In any such match the middle carbon should not have an extra substituent – that is, excluding
         the two bonds that form the pattern (to the two carbonyl carbons) it should have no extra heavy atom,
         unless it is in a ring (then one extra connectivity is tolerated).
      4. The molecule must contain a Coenzyme A adenine fragment.

    Args:
      smiles (str): SMILES string representing the molecule.

    Returns:
      bool: True if the molecule is classified as a 3-oxo-fatty acyl-CoA, False otherwise.
      str: Explanation for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject if there are any deprotonated oxygens (we want neutral CoA derivatives)
    if "[O-]" in smiles:
        return False, "Contains deprotonated oxygens; expected neutral CoA"
    
    # Check for an adenine fragment, common in Coenzyme A.
    coa_pattern = Chem.MolFromSmarts("n1cnc2ncnc(c12)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Does not contain a Coenzyme A adenine fragment"

    # Define two SMARTS patterns for the 3-oxo fatty acyl thioester fragment.
    # Pattern A forces the middle carbon to be sp3 and Pattern B allows sp2.
    frag_smarts_variants = [
        "[#6:1]-C(=O)-[#6X4:2]-C(=O)S",  # middle carbon explicitly tetrahedral
        "[#6:1]-C(=O)-[#6X3:2]-C(=O)S"   # middle carbon is trigonal (possibly in a ring)
    ]
    
    # A helper check: in the correct fragment the middle carbon (atom with tag 2 in the match)
    # should not have additional heavy-atom neighbors besides the two carbons in the pattern.
    def middle_not_over_substituted(mid_atom, patt_atom_idxs):
        # Exclude the two atoms that are part of the fragment (the two carbonyl carbons)
        exclude = set(patt_atom_idxs)
        extra_neighbors = [nbr for nbr in mid_atom.GetNeighbors() if nbr.GetIdx() not in exclude]
        # For an acetoacetyl-like fragment the middle carbon (CH2) should have no additional heavy atoms.
        # But if the carbon is in a ring, one extra connectivity (the ring closure) is acceptable.
        if mid_atom.IsInRing():
            if len(extra_neighbors) > 1:
                return False
        else:
            if len(extra_neighbors) != 0:
                return False
        return True

    valid_match_found = False
    # Try each variant of the fragment SMARTS
    for frag_smarts in frag_smarts_variants:
        frag_pattern = Chem.MolFromSmarts(frag_smarts)
        matches = mol.GetSubstructMatches(frag_pattern)
        if not matches:
            continue
        
        # For each match, further verify that the middle carbon is not over‐substituted.
        # In the substructure match,
        #   match[0] corresponds to the first carbon (R–),
        #   match[1] is the carbonyl carbon (C(=O) from the acyl chain),
        #   match[2] is the “middle” carbon (should bear only one hydrogen aside from bonds to match[1] and match[3]),
        #   match[3] is the carbonyl carbon of the thioester.
        for match in matches:
            # Safety check that there are 4 atoms in the match.
            if len(match) != 4:
                continue
            mid_atom = mol.GetAtomWithIdx(match[2])
            if not middle_not_over_substituted(mid_atom, [match[1], match[3]]):
                continue
            # Optionally, we can look for an extra carboxyl substituent on the middle carbon.
            # (For each neighbor not in the fragment, if that neighbor is a carbon that is double‐bonded to oxygen, reject.)
            extra_subst = False
            for nbr in mid_atom.GetNeighbors():
                if nbr.GetIdx() in (match[1], match[3]):
                    continue
                if nbr.GetAtomicNum() == 6:
                    for bond in nbr.GetBonds():
                        if bond.GetBondType() == Chem.BondType.DOUBLE:
                            other = bond.GetOtherAtom(nbr)
                            if other.GetAtomicNum() == 8:
                                extra_subst = True
                                break
                    if extra_subst:
                        break
            if extra_subst:
                continue
            
            # If one match clears our checks, we say we have found a valid fragment.
            valid_match_found = True
            break
        if valid_match_found:
            break

    if not valid_match_found:
        return False, "Does not contain an uncompromised 3-oxo fatty acyl thioester substructure"
    
    # Passed all checks.
    return True, "Contains a valid 3-oxo fatty acyl thioester linked to a Coenzyme A moiety"

# Example usage:
if __name__ == "__main__":
    # Example SMILES for 3-oxoadipyl-CoA:
    example_smiles = (
        "CC(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)"
        "[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)CC(=O)CCC(O)=O"
    )
    result, reason = is_3_oxo_fatty_acyl_CoA(example_smiles)
    print(result, reason)