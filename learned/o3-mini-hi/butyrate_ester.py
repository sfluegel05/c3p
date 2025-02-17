"""
Classifies: CHEBI:50477 butyrate ester
"""
"""
Classifies: butyrate ester
Definition: Any carboxylic ester where the carboxylic acid component is butyric acid.
This means that one of the ester groups in the molecule has an acyl part corresponding to n-butyrate:
CH3CH2CH2COO-
We improve upon previous attempts by using a SMARTS pattern that ensures the ester oxygen is non‐protic.
"""

from rdkit import Chem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.
    The function checks for the presence of an ester moiety wherein the acyl component is butyric acid.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule contains a butyrate ester fragment, False otherwise.
        str: A reason for the classification result.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern that represents the butyrate ester fragment.
    # Here, "[O;H0]" means an oxygen with no attached hydrogen (typical for an ester linkage) 
    # and the pattern "C(=O)CCC" covers a carbonyl attached to an ethyl chain with one extra CH2,
    # i.e. CH3CH2CH2CO- (butyrate acyl group) when matched as "O-C(=O)CCC".
    butyrate_smarts = "[O;H0]-C(=O)CCC"
    butyrate_pattern = Chem.MolFromSmarts(butyrate_smarts)
    
    # Check if the molecule has the butyrate ester substructure.
    if mol.HasSubstructMatch(butyrate_pattern):
        return True, "Found butyrate ester fragment ([O;H0]-C(=O)CCC) in the molecule."
    else:
        return False, "No butyrate ester fragment ([O;H0]-C(=O)CCC) found in the molecule."

# Example usage:
if __name__ == "__main__":
    # A couple of test SMILES strings:
    test_smiles = [
        "CCCC(=O)OCC(COC(=O)CCC)OC(=O)CCC",   # tributyrin, should be True
        "CCCC(=O)OCC",                        # ethyl butyrate, should be True
        "CCCCCCCCCCCCCCCCCC(O)C([O-])=O"        # false positive example (2-hydroxyarachidate), should be False
    ]
    for smi in test_smiles:
        result, reason = is_butyrate_ester(smi)
        print(f"SMILES: {smi}\nResult: {result} | Reason: {reason}\n")