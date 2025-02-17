"""
Classifies: CHEBI:52603 1-acyl-sn-glycero-3-phosphoserine
"""
#!/usr/bin/env python3
"""
Classifies: 1-acyl-sn-glycero-3-phosphoserine compounds
Defined as: An sn-glycerophosphoserine compound having an acyl substituent 
at the 1-hydroxy position.

Improvement rationale:
  The previous approach separately looked for a phosphoserine head group pattern 
  (P(OCC(N)C(O)=O)) and an acyl branch pattern (P(OCC(O)COC(=O)[#6])). However, some 
  false positives had extra acyl functionalities that still allowed one of these pieces 
  to match on the same phosphorus atom. To minimize this problem we now build a composite 
  SMARTS pattern that requires a phosphorus to simultaneously have:
    – a P=O (double bonded oxygen),
    – an acyl branch of the form O[CX4][CX4](O)COC(=O)[#6],
    – and a serine-containing branch of the form [O][CX4][CX4](N)C(O)=O.
  If the molecule contains a substructure matching
    P(=O)([O][CX4][CX4](N)C(O)=O)(O[CX4][CX4](O)COC(=O)[#6])
  (the order is irrelevant) then it is classified as a 1‑acyl‑sn‑glycero‑3‑phosphoserine.
"""

from rdkit import Chem

def is_1_acyl_sn_glycero_3_phosphoserine(smiles: str):
    """
    Determines whether a molecule is a 1-acyl-sn-glycero-3-phosphoserine compound.
    
    We require that the molecule contain a phosphorus center that meets all of the following:
      • It is double-bonded to an oxygen (P(=O)).
      • It bears a serine head branch attached via an oxygen; the branch should
        match [O][CX4][CX4](N)C(O)=O.
      • It bears an acyl ester branch at the 1-hydroxy position that matches
        O[CX4][CX4](O)COC(=O)[#6].
    These three groups together are encoded in the composite SMARTS pattern:
      P(=O)([O][CX4][CX4](N)C(O)=O)(O[CX4][CX4](O)COC(=O)[#6])
    We ignore chirality in the match.
    
    Args:
        smiles (str): SMILES string for the molecule.
        
    Returns:
        (bool, str): Tuple where first element is True if the compound is classified as 
                     1-acyl-sn-glycero-3-phosphoserine, else False; the second element is 
                     a reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the composite SMARTS pattern
    # Explanation:
    #  P(=O)([O][CX4][CX4](N)C(O)=O)(O[CX4][CX4](O)COC(=O)[#6])
    #   • P(=O) ensures a phosphorus atom double-bonded to an oxygen.
    #   • ([O][CX4][CX4](N)C(O)=O) captures the serine head branch.
    #   • (O[CX4][CX4](O)COC(=O)[#6]) captures the acyl branch.
    composite_smarts = "P(=O)([O][CX4][CX4](N)C(O)=O)(O[CX4][CX4](O)COC(=O)[#6])"
    composite_query = Chem.MolFromSmarts(composite_smarts)
    if composite_query is None:
        # Should never happen unless there is an error in the SMARTS syntax.
        return False, "Error in composite SMARTS pattern"
    
    # Check for serine head group separately for a better diagnostic.
    serine_smarts = "[O][CX4][CX4](N)C(O)=O"
    serine_query = Chem.MolFromSmarts(serine_smarts)
    if serine_query is None:
        return False, "Error in serine SMARTS pattern"
    
    # Check for acyl branch separately for a better diagnostic.
    acyl_smarts = "O[CX4][CX4](O)COC(=O)[#6]"
    acyl_query = Chem.MolFromSmarts(acyl_smarts)
    if acyl_query is None:
        return False, "Error in acyl branch SMARTS pattern"
    
    # First do individual checks to give specific reasons if missing
    if not mol.HasSubstructMatch(serine_query, useChirality=False):
        return False, "Phosphoserine head group not found"
    if not mol.HasSubstructMatch(acyl_query, useChirality=False):
        return False, "Acyl ester substituent at the 1-hydroxy position not found"
    
    # Now require that a phosphorus atom carries both branches in the correct overall connectivity.
    if mol.HasSubstructMatch(composite_query, useChirality=False):
        return True, "Contains phosphoserine head group with an acyl substituent at the 1-hydroxy position"
    else:
        return False, "Combined pattern of phosphoserine head group and acyl branch not found on a single phosphorus atom"


# For quick testing when running this script directly:
if __name__ == "__main__":
    test_smiles = [
        # Expected positives (from provided examples)
        "P(OCC(O)COC(=O)CCCCCCC/C=C/CCCCCCCC)(OCC(N)C(O)=O)(O)=O",   # 1-Oleoylglycerophosphoserine
        "P(OC[C@H](O)COC(=O)CCCCCCCCCCCCCCC)(OC[C@H](N)C(O)=O)(O)=O",   # PS(16:0/0:0)
        "P(OC[C@H](O)COC(=O)CCCCCCCCCCC)(OC[C@H](N)C(O)=O)(O)=O",        # PS(12:0/0:0)
        "CCCCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OC[C@H](N)C(O)=O",   # 1-stearoyl-sn-glycero-3-phosphoserine
        # Additional molecules could be tested here.
    ]
    for sm in test_smiles:
        res, reason = is_1_acyl_sn_glycero_3_phosphoserine(sm)
        print(f"SMILES: {sm}\nResult: {res}\nReason: {reason}\n")