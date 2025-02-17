"""
Classifies: CHEBI:52603 1-acyl-sn-glycero-3-phosphoserine
"""
#!/usr/bin/env python3
"""
Classifies: 1-acyl-sn-glycero-3-phosphoserine compounds
Defined as: An sn-glycerophosphoserine compound having an acyl substituent 
at the 1-hydroxy position.
Our improved algorithm checks that:
  (1) The molecule contains at least one phosphorus atom.
  (2) There is a serine branch attached to phosphorus, defined by the SMARTS
      pattern for a phosphoserine head group: P(OCC(N)C(O)=O)
  (3) There is an acyl ester substituent on phosphorus (sn-1) defined by the SMARTS
      pattern: P(OCC(O)COC(=O)[#6])
  (4) Both substructure matches share the SAME phosphorus atom (i.e. come from the same branch).
  (5) The acyl branch is found exactly once on that phosphorus atom.
If these conditions are met the molecule is classified as a 1-acyl-sn-glycero-3-phosphoserine.
If not, an appropriate reason is returned.
"""

from rdkit import Chem

def is_1_acyl_sn_glycero_3_phosphoserine(smiles: str):
    """
    Determines whether a molecule is a 1-acyl-sn-glycero-3-phosphoserine compound.
    
    The algorithm searches for a phosphorus atom that directly bears two substituents:
      - A serine head group (P(OCC(N)C(O)=O))
      - An acyl ester branch (P(OCC(O)COC(=O)[#6]))
    and ensures that the acyl branch is present exactly once.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple where the first item is True if the molecule is of this class and
                     False otherwise, and the second item is a reason explaining the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # (1) Get all phosphorus atoms in the molecule.
    phosphorus_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not phosphorus_atoms:
        return False, "No phosphorus atom found"
    
    # (2) Define SMARTS patterns that are anchored on phosphorus.
    # The serine branch is expected to be: P(OCC(N)C(O)=O)
    serine_smarts = "P(OCC(N)C(O)=O)"
    serine_query = Chem.MolFromSmarts(serine_smarts)
    if serine_query is None:
        return False, "Error in serine SMARTS pattern"

    # The acyl branch is expected to be: P(OCC(O)COC(=O)[#6])
    acyl_smarts = "P(OCC(O)COC(=O)[#6])"
    acyl_query = Chem.MolFromSmarts(acyl_smarts)
    if acyl_query is None:
        return False, "Error in acyl SMARTS pattern"

    # (3) Find all substructure matches for each pattern (they are tuples with the phosphorus atom first).
    serine_matches = mol.GetSubstructMatches(serine_query, useChirality=False)
    acyl_matches = mol.GetSubstructMatches(acyl_query, useChirality=False)

    if not serine_matches:
        return False, "Phosphoserine head group not found"
    if not acyl_matches:
        return False, "Acyl ester substituent at the 1-hydroxy position not found"

    # (4) We now check whether there is at least one phosphorus atom for which
    # both the serine and acyl branch are attached.
    # In our SMARTS both patterns have the P at index 0.
    serine_P = {match[0] for match in serine_matches}
    acyl_P = {match[0] for match in acyl_matches}
    common_P = serine_P.intersection(acyl_P)
    if not common_P:
        return False, "No single phosphorus atom has both a serine head group and an acyl branch"

    # (5) For each phosphorus in common_P, enforce that the acyl branch occurs exactly once.
    for p_idx in common_P:
        acyl_count = sum(1 for match in acyl_matches if match[0] == p_idx)
        serine_count = sum(1 for match in serine_matches if match[0] == p_idx)
        if serine_count < 1:
            continue  # Should not happen because of the intersection test.
        if acyl_count == 1:
            return True, "Contains phosphoserine head group with an acyl substituent at the 1-hydroxy position"
        else:
            return False, f"Found {acyl_count} acyl ester substituents on the phosphorus atom; expected exactly 1"

    return False, "No phosphorus substituent with valid serine and acyl branches found"


# Quick testing when running this script directly
if __name__ == "__main__":
    test_smiles = [
        # Expected positives
        "P(OCC(O)COC(=O)CCCCCCC/C=C/CCCCCCCC)(OCC(N)C(O)=O)(O)=O",  # 1-Oleoylglycerophosphoserine
        "P(OC[C@H](O)COC(=O)CCCCCCCCCCCCCCC)(OC[C@H](N)C(O)=O)(O)=O",  # PS(16:0/0:0)
        "P(OC[C@H](O)COC(=O)CCCCCCCCCCC)(OC[C@H](N)C(O)=O)(O)=O",       # PS(12:0/0:0)
        "CCCCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OC[C@H](N)C(O)=O",  # 1-stearoyl-sn-glycero-3-phosphoserine
        # An example false positive from previous attempts might be added here for testing.
    ]
    for sm in test_smiles:
        res, reason = is_1_acyl_sn_glycero_3_phosphoserine(sm)
        print(f"SMILES: {sm}\nResult: {res}\nReason: {reason}\n")