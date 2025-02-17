"""
Classifies: CHEBI:23003 carbamate ester
"""
"""
Classifies: Carbamate ester 
Defined as: Any ester of carbamic acid or its N-substituted derivatives.
Functional group: R–O–C(=O)–N
"""

from rdkit import Chem

def is_carbamate_ester(smiles: str):
    """
    Determines if a molecule is a carbamate ester based on its SMILES string.
    A carbamate ester is defined as an ester of carbamic acid (or an N-substituted derivative)
    having the functional group R–O–C(=O)–N.
    
    This function uses a SMARTS pattern that explicitly represents the core
    functional group as four atoms:
        - an arbitrary atom (R)
        - an oxygen atom (O)
        - a carbonyl carbon (C) with an explicit double bond to oxygen (=[O])
        - a nitrogen atom (N)
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        (bool, str): Tuple where the boolean is True if the molecule is classified as a carbamate ester,
                     and the string provides a reason. Otherwise, returns False and an error message.
    """
    
    # Convert SMILES to an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for the carbamate ester functional group: R-O-C(=O)-N.
    # The pattern uses explicit atom mapping:
    # [*:1] matches any atom (the R group),
    # [O:2] matches the oxygen connecting the R to the carbonyl,
    # [C:3](=[O]) matches the carbonyl carbon (with a double bond to an oxygen),
    # [N:4] matches the nitrogen.
    smarts = "[*:1]-[O:2]-[C:3](=[O])-[N:4]"
    pattern = Chem.MolFromSmarts(smarts)
    if pattern is None:
        return False, "Error creating SMARTS pattern"
    
    # Find all substructure matches of the pattern in the molecule
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "Carbamate ester pattern not found in the molecule"
    
    # If at least one match is found, we assume the motif is present.
    # (The SMARTS explicitly requires the double bond to oxygen, which ensures that
    # the connectivity is as expected.)
    for match in matches:
        if len(match) == 4:
            return True, "Contains carbamate ester functional group (R-O-C(=O)-NR')"
    
    # If no match passed the check, report a connectivity issue.
    return False, "Carbamate ester pattern found but did not meet connectivity requirements"

# Example usage (for testing; these lines can be commented out or removed in deployment)
if __name__ == "__main__":
    test_smiles = [
        "COc1ccccc1OCC(O)COC(N)=O",  # 2-hydroxy-3-(2-methoxyphenoxy)propyl carbamate (should be True)
        "CN(C)C(=O)Oc1cc(OC(=O)N(C)C)cc(c1)C(O)CNC(C)(C)C",  # bambuterol (should be True)
        "S(SCCNC(=O)CC[C@H](N)C(O)=O)CCNC=1C(=O)C=2[C@H]([C@]3(OC)N(C[C@@H]4N[C@H]34)C2C(=O)C1C)COC(=O)N",  # mitomycin derivative (should be True)
        "C(=O)(OCC1=CN=CC=C1)NCC2=CC=C(C=C2)C(=O)NC3=CC=C(C=C3N)F",  # pyridin-3-ylmethyl carbamate (should be True)
        "CCOC(N)=O",  # urethane (should be True)
        "CNC(=O)Oc1cccc2ccccc12",  # carbaryl (should be True)
    ]
    for smi in test_smiles:
        res, reason = is_carbamate_ester(smi)
        print(f"SMILES: {smi}\nResult: {res}\nReason: {reason}\n")