"""
Classifies: CHEBI:23003 carbamate ester
"""
"""
Classifies: Carbamate ester 
Defined as: Any ester of carbamic acid or its N-substituted derivatives.
The functional group is characterized as R–O–C(=O)–NR.
"""

from rdkit import Chem

def is_carbamate_ester(smiles: str):
    """
    Determines if a molecule is a carbamate ester based on its SMILES string.
    A carbamate ester is defined as an ester of carbamic acid (or an N-substituted derivative)
    having the functional group R–O–C(=O)–NR.

    This implementation does two steps:
      1. It uses a SMARTS pattern "[*]-O-C(=O)-[NX2,NX3]" to search for substructures matching 
         an arbitrary R attached to an ester oxygen, a carbonyl function and a nitrogen.
      2. For each match, it relaxes the connectivity check: the carbonyl carbon (C) must have
         (i) the ester oxygen (O) and (ii) the nitrogen (N) among its neighbors, 
         and (iii) at least one oxygen in a double bond (the carbonyl oxygen) not explicitly 
         captured in the SMARTS. In this way, additional substituents or ring closures are allowed.
         
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        (bool, str): Tuple with boolean True if the molecule is classified as a carbamate ester,
                     and a corresponding reason. Otherwise, False and the reason why it failed.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS to capture the carbamate ester motif: R-O-C(=O)-N
    smartsmarts = "[*:1]-O:2-C(=O)-[NX2,NX3:3]"
    pattern = Chem.MolFromSmarts(smartsmarts)
    if pattern is None:
        return False, "Error creating SMARTS pattern"
    
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "Carbamate ester pattern not found in the molecule"
    
    # Loop over all substructure matches to check for proper connectivity.
    for match in matches:
        # The SMARTS was defined so that match = (r_idx, o_idx, c_idx, n_idx)
        if len(match) != 4:
            continue  # skip if not the expected number of atoms
        
        r_idx, o_idx, c_idx, n_idx = match
        c_atom = mol.GetAtomWithIdx(c_idx)
        # Get the neighbors of the carbonyl carbon
        neighbor_idxs = [nbr.GetIdx() for nbr in c_atom.GetNeighbors()]
        
        # Ensure that the ester oxygen and the nitrogen are in the neighbors.
        if o_idx not in neighbor_idxs or n_idx not in neighbor_idxs:
            continue
        
        # Now, check that the carbonyl carbon is double-bonded to at least one oxygen
        # (which has not been explicitly matched as the ester oxygen).
        has_carbonyl_oxygen = False
        for nbr in c_atom.GetNeighbors():
            if nbr.GetIdx() == o_idx:
                continue  # skip the ester oxygen from the match
            bond = mol.GetBondBetweenAtoms(c_idx, nbr.GetIdx())
            if bond is not None and bond.GetBondTypeAsDouble() == 2 and nbr.GetAtomicNum() == 8:
                has_carbonyl_oxygen = True
                break
        
        if has_carbonyl_oxygen:
            return True, "Contains carbamate ester functional group (R-O-C(=O)-NR')"
    
    # If no match passes the connectivity check, then the molecule is not classified as a carbamate ester.
    return False, "Carbamate ester pattern found but did not meet connectivity requirements"

# Example usage (for testing; can be removed or commented out):
if __name__ == "__main__":
    test_smiles = [
        "COc1ccccc1OCC(O)COC(N)=O",  # 2-hydroxy-3-(2-methoxyphenoxy)propyl carbamate (should be True)
        "CN(C)C(=O)Oc1cc(OC(=O)N(C)C)cc(c1)C(O)CNC(C)(C)C",  # bambuterol (should be True)
        "S(SCCNC(=O)CC[C@H](N)C(O)=O)CCNC=1C(=O)C=2[C@H]([C@]3(OC)N(C[C@@H]4N[C@H]34)C2C(=O)C1C)COC(=O)N",  # mitomycin derivative (should be True)
        "C(=O)(OCC1=CN=CC=C1)NCC2=CC=C(C=C2)C(=O)NC3=CC=C(C=C3N)F",  # pyridin-3-ylmethyl ...carbamate (should be True)
        "CCOC(N)=O",  # urethane (should be True)
        "CNC(=O)Oc1cccc2ccccc12",  # carbaryl (should be True)
    ]
    for smi in test_smiles:
        res, reason = is_carbamate_ester(smi)
        print(f"SMILES: {smi}\nResult: {res}\nReason: {reason}\n")