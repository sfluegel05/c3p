"""
Classifies: CHEBI:23003 carbamate ester
"""
"""
Classifies: Carbamate ester 
Defined as: Any ester of carbamic acid or its N‐substituted derivatives.
The characteristic functional group is R–O–C(=O)–NR′.
"""

from rdkit import Chem

def is_carbamate_ester(smiles: str):
    """
    Determines if a molecule is a carbamate ester based on its SMILES string.
    A carbamate ester has the functional group R–O–C(=O)–NR, meaning an oxygen connected
    to a carbonyl carbon that is in turn bound to a nitrogen (which may be aromatic or aliphatic).

    This function uses a more refined SMARTS pattern designed to ensure:
      - The oxygen comes from an R-group (thus the pattern "[#6]-O"),
      - The carbon is part of a carbonyl (C(=O)),
      - And the bonded nitrogen can be either aliphatic (N) or aromatic (n).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a carbamate ester, False otherwise.
        str: Reason for the classification.
    """

    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a refined SMARTS pattern for a carbamate ester:
    # "[#6]-O-C(=O)-[$([NX3]),$([nX3])]"
    #   [#6]         : any carbon, ensuring the oxygen comes from an R-group.
    #   O            : the ester oxygen.
    #   C(=O)       : the carbonyl group.
    #   [$([NX3]),$([nX3])] : a nitrogen that is either aliphatic (N) or aromatic (n).
    carbamate_pattern = Chem.MolFromSmarts("[#6]-O-C(=O)-[$([NX3]),$([nX3])]")
    if carbamate_pattern is None:
        return False, "Error in SMARTS pattern creation"
    
    # Check if the molecule contains at least one match of the carbamate ester substructure.
    matches = mol.GetSubstructMatches(carbamate_pattern)
    if matches:
        return True, "Contains carbamate ester functional group (R-O-C(=O)-NR')"
    else:
        return False, "Carbamate ester pattern not found in the molecule"

# Example usage (for testing purposes; these lines can be commented out or removed):
if __name__ == "__main__":
    # Test examples (both positives and known negatives)
    test_smiles = [
        "COc1ccccc1OCC(O)COC(N)=O",       # 2-hydroxy-3-(2-methoxyphenoxy)propyl carbamate (positive)
        "CN(C)C(=O)Oc1cc(OC(=O)N(C)C)cc(c1)C(O)CNC(C)(C)C",  # bambuterol (positive)
        "CCOC(N)=O",                     # urethane (positive)
        "CNC(=O)Oc1cccc2ccccc12",         # carbaryl (positive)
        "C(CN(C)C)(=O)N1C=2C(CCC3=C1C=CC=C3)=CC=C(C2)NC(OCC)=O",  # tiracizine (false positive candidate)
        "CCOC(=O)n1ccn(C)c1=S",           # carbimazole (false negative candidate)
    ]
    for smi in test_smiles:
        result, reason = is_carbamate_ester(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")