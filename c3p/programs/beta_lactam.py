"""
Classifies: CHEBI:35627 beta-lactam
"""
"""
Classifies: beta-lactam
A beta-lactam is defined as a cyclic amide in which the amide (nitrogen–carbonyl) bond is
contained within a 4-membered ring.
This implementation:
  - Parses the SMILES, selects the largest fragment and adds explicit hydrogens.
  - Uses a SMARTS query to target the core beta-lactam motif:
        [NX3;R4][CX3;R4](=O)[CX3;R4][CX3;R4]
    This pattern requires a 4-membered ring (R4) that contains exactly one nitrogen (NX3)
    and three carbons (CX3), with one of the carbon atoms bearing an external carbonyl oxygen.
  - For each substructure match the code further verifies that:
      (a) the oxygen double-bond is outside the ring,
      (b) the bond between the nitrogen and the carbonyl carbon is a single bond.
If such a substructure is found, the molecule is classified as a beta-lactam.
"""

from rdkit import Chem

def is_beta_lactam(smiles: str):
    """
    Determines if a molecule is a beta-lactam based on its SMILES string.

    The function uses a SMARTS query for the core beta-lactam ring, defined as a 
    4-membered ring consisting of one nitrogen and three carbons where one carbon (the carbonyl)
    is double bonded to an oxygen outside of the ring and is directly bonded to the ring nitrogen 
    through a single bond (the amide bond).

    Args:
        smiles (str): SMILES string representing the molecule.
    
    Returns:
        bool: True if the molecule is classified as a beta-lactam, False otherwise.
        str: A reason supporting the classification decision.
    """
    # Parse the molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all fragments and choose the largest to avoid solvents/salts.
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return False, "No fragments could be parsed from the SMILES"
    main_frag = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    
    # Add explicit hydrogens to ensure proper bond perception.
    main_frag = Chem.AddHs(main_frag)
    
    # Define a SMARTS query for a candidate beta-lactam ring:
    # - The "R4" flag forces ring membership in a 4-membered ring.
    # - The pattern requires a ring-bound trivalent nitrogen,
    #   followed by a ring-bound carbon that is double-bonded to an oxygen,
    #   followed by two ring-bound carbons.
    beta_lactam_smarts = "[NX3;R4][CX3;R4](=O)[CX3;R4][CX3;R4]"
    query = Chem.MolFromSmarts(beta_lactam_smarts)
    if query is None:
        return False, "Error parsing SMARTS for beta-lactam"
    
    # Find all substructure matches in the main fragment.
    matches = main_frag.GetSubstructMatches(query)
    if not matches:
        return False, "No beta-lactam ring found (4-membered ring with the required amide bond was not detected)"
    
    # Iterate over each match and perform additional checks.
    for match in matches:
        # match is a tuple of atom indices corresponding to the SMARTS ordering:
        # [N_idx, C_carbonyl_idx, C_idx, C_idx]
        N_idx, C_carbonyl_idx, _, _ = match
        
        # Fetch the bond between the nitrogen and the carbonyl carbon.
        bond = main_frag.GetBondBetweenAtoms(N_idx, C_carbonyl_idx)
        if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
            # If the bond linking the N and the carbonyl carbon is not single, skip this match.
            continue
        
        # Now verify that the oxygen (double-bonded to the carbonyl carbon) is not part of the ring.
        C_carbonyl = main_frag.GetAtomWithIdx(C_carbonyl_idx)
        oxygen_found = False
        for nbr in C_carbonyl.GetNeighbors():
            # Skip if neighbor is part of the match (i.e. in the beta-lactam ring)
            if nbr.GetIdx() in match:
                continue
            if nbr.GetAtomicNum() == 8:
                ext_bond = main_frag.GetBondBetweenAtoms(C_carbonyl_idx, nbr.GetIdx())
                if ext_bond is not None and ext_bond.GetBondType() == Chem.BondType.DOUBLE:
                    oxygen_found = True
                    break
        if not oxygen_found:
            continue  # no valid external carbonyl oxygen found in this match
        
        # If all criteria are met for this substructure match, then we have a beta-lactam.
        return True, ("Beta-lactam ring detected: found a 4-membered ring with one nitrogen and three carbons, "
                      "where one carbon bears an external double-bonded oxygen and is directly bonded to the nitrogen "
                      "by a single bond representing the amide linkage")
    
    # If no substructure match completely fulfills all criteria, return negative.
    return False, "No beta-lactam ring found (none of the candidate rings fulfilled all criteria)"

# Example usage:
if __name__ == "__main__":
    # A few example SMILES strings for testing.
    test_smiles = [
        "O=C1CCN1",  # azetidin-2-one, simplest beta-lactam
        "[H][C@]12SCC(C)=C(N1C(=O)[C@H]2N)C(O)=O",  # 7beta-aminodeacetoxycephalosporanic acid
        "C[C@H]1[C@@H](C(=O)N1S(=O)(=O)O)NC(=O)C(=NOC(C)(C)C(=O)O)C2=CSC(=N2)N"  # complex beta-lactam derivative
    ]
    
    for smi in test_smiles:
        result, reason = is_beta_lactam(smi)
        print("SMILES:", smi)
        print("Result:", result)
        print("Reason:", reason)
        print("------")