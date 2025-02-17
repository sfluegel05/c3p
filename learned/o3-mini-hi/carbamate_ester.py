"""
Classifies: CHEBI:23003 carbamate ester
"""
"""
Classifies: Carbamate ester 
Defined as: Any ester of carbamic acid or its N‐substituted derivatives.
The functional group is characterized as R–O–C(=O)–NR′.
"""

from rdkit import Chem

def is_carbamate_ester(smiles: str):
    """
    Determines if a molecule is a carbamate ester based on its SMILES string.
    A carbamate ester has the functional group R–O–C(=O)–NR, indicating an 
    ester oxygen attached to a carbonyl carbon that is in turn bound to a nitrogen.

    This implementation uses a two‐step approach:
      1. A SMARTS pattern search for any substructure of the form
         “[*]-O-C(=O)-[NX2,NX3]” (i.e. any heavy atom “R” attached to O, then C(=O),
         then a nitrogen that can be sp2 or sp3).
      2. A post‐match check on the carbonyl carbon’s local bonding – it should have
         exactly three neighbors: (i) the ester oxygen from the substructure,
         (ii) a double‐bonded oxygen (the carbonyl oxygen) and (iii) the nitrogen.
         This helps avoid false positives from “accidental” matches.
         
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a carbamate ester, False otherwise.
        str: Reason for the classification.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS: any heavy atom (*) attached to an O, then a carbonyl C(=O), then a nitrogen.
    # Note: We allow the nitrogen to be either sp2 or sp3 ([NX2,NX3]).
    pattern = Chem.MolFromSmarts("[*]-O-C(=O)-[NX2,NX3]")
    if pattern is None:
        return False, "Error in SMARTS pattern creation"
    
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "Carbamate ester pattern not found in the molecule"
    
    # Post-match filtering: for each match, check that the C in C(=O) has exactly 3 neighbors:
    #   (a) the ester oxygen (match[0])
    #   (b) the nitrogen (match[2])
    #   (c) a double-bonded oxygen (carbonyl oxygen) that is NOT the one from the match.
    for match in matches:
        # match returns a tuple of atom indices [r, O, C, N] as we defined the pattern order:
        # The SMARTS "[*]-O-C(=O)-[NX2,NX3]" gives:
        #   idx0: the atom attached to O (R, any heavy atom)
        #   idx1: the ester oxygen
        #   idx2: the carbonyl carbon
        #   idx3: the nitrogen
        if len(match) != 4:
            # Our pattern defined 4 atoms; if not, skip
            continue
        
        r_idx, o_idx, c_idx, n_idx = match
        c_atom = mol.GetAtomWithIdx(c_idx)
        
        # The carbonyl carbon must have exactly 3 neighbors.
        if c_atom.GetDegree() != 3:
            continue
        
        neighbors = [nbr.GetIdx() for nbr in c_atom.GetNeighbors()]
        # It must be bonded to both the ester oxygen and the nitrogen from the match.
        if o_idx not in neighbors or n_idx not in neighbors:
            continue
        
        # Identify the carbonyl oxygen: it should be a neighbor (different from the ester oxygen)
        carbonyl_os = []
        for nbr in c_atom.GetNeighbors():
            if nbr.GetIdx() == o_idx:
                continue  # skip the ester oxygen
            bond = mol.GetBondBetweenAtoms(c_idx, nbr.GetIdx())
            # Check if this bond is a double bond and neighbor is oxygen.
            if bond is not None and bond.GetBondTypeAsDouble() == 2 and nbr.GetAtomicNum() == 8:
                carbonyl_os.append(nbr.GetIdx())
        if len(carbonyl_os) != 1:
            # The carbonyl carbon must have exactly one double-bonded oxygen (carbonyl oxygen)
            continue
        
        # If all conditions are met, then we consider the molecule to contain a valid carbamate ester.
        return True, "Contains carbamate ester functional group (R-O-C(=O)-NR')"
    
    # If no valid substructure passes the filtering, then it is not classified as a carbamate ester.
    return False, "Carbamate ester pattern found but did not meet connectivity requirements"

# Example usage (for testing purposes; can be commented out or removed):
if __name__ == "__main__":
    test_smiles = [
        "COc1ccccc1OCC(O)COC(N)=O",        # 2-hydroxy-3-(2-methoxyphenoxy)propyl carbamate (true positive)
        "CN(C)C(=O)Oc1cc(OC(=O)N(C)C)cc(c1)C(O)CNC(C)(C)C",  # bambuterol (true positive)
        "CCOC(N)=O",                      # urethane (true positive)
        "CNC(=O)Oc1cccc2ccccc12",          # carbaryl (true positive)
        "CNC(=O)ON=C(SC)C(=O)N(C)C",        # oxamyl (previously false negative candidate)
        "CNC(=O)O\\N=C(\\C)C(C)SC"          # butocarboxim (previously false negative candidate)
    ]
    for smi in test_smiles:
        res, reason = is_carbamate_ester(smi)
        print(f"SMILES: {smi}\nResult: {res}\nReason: {reason}\n")