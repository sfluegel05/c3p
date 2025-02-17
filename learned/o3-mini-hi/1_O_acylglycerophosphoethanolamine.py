"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
"""
Classifies: 1-O-acylglycerophosphoethanolamine
Definition: A glycerophosphoethanolamine having an unspecified O-acyl substituent 
at the 1-position of the glycerol fragment.

This improved heuristic checks for the full combined motif:
  (1) A glycerol backbone bound at one end to a phosphate that is in turn attached to an ethanolamine headgroup.
  (2) The glycerol backbone is required to have an ester linkage (i.e. an acyl group attached via a C(=O)O– bond)
      on the sn-1 position.
      
The SMARTS below (ignoring stereochemistry) 
   O[CH2][CH](O)[CH2]O[P](=O)(O)OCCN
matches a glycerophosphoethanolamine (GPE) motif in which the glycerol is represented as three carbons 
(implicitly CH2-CH(OH)-CH2) and the phosphate is attached to the terminal CH2.
Once this GPE motif is found, we identify the oxygen immediately preceding the sn-1 CH2 
(i.e. the first “O” in the SMARTS) and require that it be part of an ester bond:
that is, it must be connected to a carbon (the acyl carbon) that is doubly bonded to an oxygen.
"""

from rdkit import Chem

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine based on its SMILES string.
    
    Heuristic:
      1. The molecule must contain a glycerophosphoethanolamine (GPE) headgroup.
         We require the substructure matching a glycerol fragment attached to a phosphate and ethanolamine:
           O[CH2][CH](O)[CH2]O[P](=O)(O)OCCN
      2. In that fragment, the oxygen preceding the sn-1 CH2 must come from an ester bond.
         That is, that oxygen (call it ester-O) should be bonded to an acyl carbon (C)
         which in turn is double-bonded to an oxygen (a carbonyl).
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if the molecule matches this 1-O-acylglycerophosphoethanolamine motif, else False.
      str: Explanation of the decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for the glycerophosphoethanolamine backbone.
    # The pattern (ignoring stereochemistry) is:
    #   A terminal oxygen (of the ester linkage),
    #   followed by a CH2 (sn-1 carbon),
    #   then a CH with an OH (sn-2),
    #   then a CH2,
    #   then an oxygen connected to phosphorus which is double-bonded to an O and has O attached to CCN.
    gpe_smarts = "O[CH2][CH](O)[CH2]O[P](=O)(O)OCCN"
    gpe_pattern = Chem.MolFromSmarts(gpe_smarts)
    if gpe_pattern is None:
        return False, "Error in GPE SMARTS pattern"
    
    matches = mol.GetSubstructMatches(gpe_pattern)
    if not matches:
        return False, "Glycerophosphoethanolamine headgroup not identified"

    # Now, for each GPE match, check the acyl substitution.
    # According to our SMARTS, the atoms in the match are (by order):
    #   idx0: The oxygen that should come from the ester linkage (ester_O).
    #   idx1: The sn-1 CH2 carbon.
    #   idx2: The sn-2 CH bearing the free OH.
    #   idx3: The sn-3 CH2 (attached to phosphate).
    #   idx4: The oxygen linking glycerol to phosphorus.
    #   idx5: The phosphorus atom.
    #   idx6, idx7, idx8: The remainder of the ethanolamine headgroup (O, C, C, N... as applicable).
    # We then verify that the ester_O (atom idx0) is connected to an acyl carbon that is in a carbonyl.
    for match in matches:
        # Get the atom index for the ester oxygen from our pattern:
        ester_oxygen_idx = match[0]
        sn1_idx = match[1]
        ester_oxygen = mol.GetAtomWithIdx(ester_oxygen_idx)
        
        # Look at neighbors of the ester oxygen.
        # One neighbor should be the sn-1 carbon (we already know that) and the other should be from the acyl chain.
        acyl_carbon = None
        for nbr in ester_oxygen.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx == sn1_idx:
                continue  # skip the glycerol part
            # Otherwise assume this neighbor is from the acyl chain:
            acyl_carbon = nbr
            break
        if acyl_carbon is None:
            continue  # no acyl chain detected in this match; check next match

        # Check that the acyl_carbon is indeed in a carbonyl environment,
        # i.e. that it has at least one double bond to an oxygen.
        has_carbonyl = False
        for bond in acyl_carbon.GetBonds():
            if bond.GetBondTypeAsDouble() == 2:  # double bond?
                other_atom = bond.GetOtherAtom(acyl_carbon)
                if other_atom.GetAtomicNum() == 8:  # oxygen
                    has_carbonyl = True
                    break
        if has_carbonyl:
            return True, "Molecule contains a phosphoethanolamine headgroup with 1-O-acyl substitution on glycerol"
    
    return False, "Glycerophosphoethanolamine headgroup found but acyl substitution at sn-1 not clearly identified"


# Example usage (uncomment to try):
# smiles_example = "CCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OCCN"  # 1-heptadecanoyl-sn-glycero-3-phosphoethanolamine
# result, reason = is_1_O_acylglycerophosphoethanolamine(smiles_example)
# print(result, reason)