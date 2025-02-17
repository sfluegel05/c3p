"""
Classifies: CHEBI:50477 butyrate ester
"""
"""
Classifies: butyrate ester
Definition: Any carboxylic ester where the carboxylic acid component is butyric acid.
That is, the acyl part is butyrate, CH3CH2CH2CO-, so the ester fragment is -O-C(=O)[CH2][CH2][CH3].
We refine the SMARTS pattern by ensuring:
  • The ester oxygen is non‐protic and not anionic.
  • The butyrate acyl group is exactly three carbons long ending in a terminal methyl.
"""

from rdkit import Chem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.
    We require that the molecule contains an ester moiety in which the acyl portion is butyric acid,
    that is, the substructure should be –O–C(=O)–CH2–CH2–CH3 (with the terminal methyl having only 1 neighbor).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule contains a butyrate ester fragment, False otherwise.
        str: A reason explaining the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Define an improved SMARTS pattern for the butyrate ester fragment.
    # Explanation:
    #   [O;H0;!$([O-])]  Match an oxygen atom with no attached hydrogen (H0) and it is not anionic.
    #   -C(=O)          The oxygen is bonded to a carbonyl carbon.
    #   [CH2]           Next, a methylene group.
    #   [CH2]           Followed by another methylene group.
    #   [CH3;D1]        Finally, a methyl group that has a degree of 1 (i.e. terminal).
    butyrate_smarts = "[O;H0;!$([O-])]-C(=O)[CH2][CH2][CH3;D1]"
    
    butyrate_pattern = Chem.MolFromSmarts(butyrate_smarts)
    if butyrate_pattern is None:
        return False, "Invalid SMARTS pattern."
    
    # Check if the molecule has a substructure match for the butyrate ester fragment.
    if mol.HasSubstructMatch(butyrate_pattern):
        return True, "Found butyrate ester fragment ([O;H0;!$([O-])]-C(=O)[CH2][CH2][CH3;D1]) in the molecule."
    else:
        return False, "No butyrate ester fragment ([O;H0;!$([O-])]-C(=O)[CH2][CH2][CH3;D1]) found in the molecule."

# Example usage: (for testing purposes)
if __name__ == "__main__":
    test_smiles = [
        "CCCC(=O)OCC(COC(=O)CCC)OC(=O)CCC",  # tributyrin, should be True
        "CCCC(=O)OCC",                       # ethyl butyrate, should be True
        "CCCCCCCCCCCCCCCCCC(O)C([O-])=O"       # false positive example (2-hydroxyarachidate), should be False
    ]
    for smi in test_smiles:
        result, reason = is_butyrate_ester(smi)
        print(f"SMILES: {smi}\nResult: {result} | Reason: {reason}\n")