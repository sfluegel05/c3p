"""
Classifies: CHEBI:23003 carbamate ester
"""
"""
Classifies: Carbamate ester
Defined as: Any ester of carbamic acid or its N-substituted derivatives.
Functional group: R–O–C(=O)–N

This revised implementation uses a SMARTS pattern without enforcing an exact match length.
If the substructure is found anywhere in the molecule, we consider the molecule to be a carbamate ester.
"""

from rdkit import Chem

def is_carbamate_ester(smiles: str):
    """
    Determines if a molecule is a carbamate ester based on its SMILES string.
    A carbamate ester is defined as an ester of carbamic acid (or an N-substituted derivative)
    having the functional group R–O–C(=O)–N.

    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        (bool, str): A tuple where the boolean is True if the molecule contains a carbamate ester
                     functional group and the string explains the reason. Otherwise, returns False.
    """
    
    # Convert SMILES to an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for the carbamate ester functional group: R-O-C(=O)-N.
    # Here the pattern is set as a simple anchor: any atom (R) connected to an oxygen,
    # which is connected to a carbonyl carbon (C with a double bond to O) and then to a nitrogen.
    # We do not explicitly enforce that the substructure match has length 4 in order to allow
    # for extensions due to substituents.
    smarts = "[*:1]-O-C(=O)-N"
    pattern = Chem.MolFromSmarts(smarts)
    if pattern is None:
        return False, "Error creating SMARTS pattern"
    
    # Find all substructure matches of the pattern in the molecule
    matches = mol.GetSubstructMatches(pattern)
    if matches:
        return True, "Contains carbamate ester functional group (R-O-C(=O)-N)"
    
    return False, "Carbamate ester pattern not found in the molecule"

# The following example usage lines are for testing purposes 
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