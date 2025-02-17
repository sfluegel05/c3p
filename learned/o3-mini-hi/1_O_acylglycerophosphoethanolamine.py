"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
"""
Classifies: 1-O-acylglycerophosphoethanolamine
Definition: A glycerophosphoethanolamine having an unspecified O-acyl substituent 
at the 1-position of the glycerol fragment.
Heuristic:
  (1) The molecule must contain a glycerophosphoethanolamine (GPE) headgroup.
      We demand a substructure matching a glycerol fragment esterified at sn-1 (via oxygen)
      and linked to a phosphate that carries an ethanolamine (must match exactly “OCCN”).
      The SMARTS used is:
           O[CH2][CH](O)[CH2]O[P](=O)(O)OCCN
      which assumes the glycerol is represented as: O-CH2-CH(OH)-CH2-O-
      followed by the phosphate: P(=O)(O) attached to an ethanolamine headgroup: OCCN.
  (2) In that match the first atom (an oxygen) should be the one bridging an acyl group.
      We then require that this oxygen (the ester-O) be connected (besides the sn-1 CH2)
      to a carbon that is part of a carbonyl group (i.e. has a double bond to oxygen).
  (3) We also require that the ester oxygen have exactly two neighbors
      (one to the glycerol and the other to the acyl carbon)
      so as to avoid spurious matches.
If these are satisfied, we return True.
"""

from rdkit import Chem

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine based on its SMILES string.
    
    Args:
      smiles (str): SMILES string of the molecule.

    Returns:
      bool: True if the molecule is classified as a 1-O-acylglycerophosphoethanolamine, else False.
      str: Explanation of the decision.
    """
    # Parse SMILES into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for the glycerophosphoethanolamine backbone.
    # The pattern (ignoring stereochemistry) is:
    #   O[CH2][CH](O)[CH2]O[P](=O)(O)OCCN
    # Here the very first "O" is the putative ester oxygen (ester_O) that should be part of an
    # acyl linkage at the glycerol sn-1 position.
    gpe_smarts = "O[CH2][CH](O)[CH2]O[P](=O)(O)OCCN"
    gpe_pattern = Chem.MolFromSmarts(gpe_smarts)
    if gpe_pattern is None:
        return False, "Error in constructing SMARTS for GPE pattern"
    
    matches = mol.GetSubstructMatches(gpe_pattern)
    if not matches:
        return False, "Glycerophosphoethanolamine headgroup not identified"

    # For each match of the GPE headgroup, verify the acyl substitution on the sn-1 oxygen.
    # Based on our SMARTS the match indices are assumed (in order) as:
    #   idx0: ester oxygen (expected to be connected to an acyl carbon as well as the sn-1 CH2)
    #   idx1: sn-1 CH2 (glycerol)
    #   idx2: sn-2 CH with free OH
    #   idx3: sn-3 CH2 (attached to the phosphate)
    #   idx4: oxygen that links to phosphorus
    #   idx5: phosphorus
    #   idx6, idx7, idx8: atoms making up the ethanolamine (O, C, C, N) – where we demand "OCCN"
    for match in matches:
        ester_oxygen_idx = match[0]
        sn1_idx = match[1]
        ester_oxygen = mol.GetAtomWithIdx(ester_oxygen_idx)

        # Enforce that the candidate ester oxygen must have exactly two neighbors.
        if ester_oxygen.GetDegree() != 2:
            continue

        # Look among the neighbors: one must be the sn-1 carbon (already identified) and the other should be
        # the acyl chain (expected to be a carbon that bears a carbonyl group).
        acyl_carbon = None
        for nbr in ester_oxygen.GetNeighbors():
            if nbr.GetIdx() == sn1_idx:
                continue  # Skip the glycerol connection.
            acyl_carbon = nbr
            break
        if acyl_carbon is None:
            continue  # No candidate acyl carbon found in this match.

        # Check that the acyl carbon is indeed part of a carbonyl: at least one of its bonds
        # must be a double bond to an oxygen.
        has_carbonyl = False
        # Iterate over bonds of the acyl carbon:
        for bond in acyl_carbon.GetBonds():
            # bond.GetBondTypeAsDouble() returns 2.0 if it is a double bond.
            if bond.GetBondTypeAsDouble() == 2:
                other = bond.GetOtherAtom(acyl_carbon)
                if other.GetAtomicNum() == 8:  # oxygen
                    has_carbonyl = True
                    break
        if has_carbonyl:
            return True, ("Molecule contains a phosphoethanolamine headgroup with 1-O-acyl substitution on glycerol, "
                          "with the ester oxygen showing proper connectivity to a carbonyl-bearing acyl chain.")
    return False, ("Glycerophosphoethanolamine headgroup found but acyl substitution at sn-1 not clearly identified, "
                   "or the ester oxygen does not display the expected bonding environment.")
    
# Example usage:
# Uncomment the following lines to test with one of the provided SMILES:
# smiles_example = "CCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OCCN"  # 1-heptadecanoyl-sn-glycero-3-phosphoethanolamine
# result, reason = is_1_O_acylglycerophosphoethanolamine(smiles_example)
# print(result, reason)