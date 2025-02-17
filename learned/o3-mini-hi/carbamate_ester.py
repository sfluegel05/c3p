"""
Classifies: CHEBI:23003 carbamate ester
"""
"""
Classifies: Carbamate ester 
Defined as: Any ester of carbamic acid or its N‐substituted derivatives.
The characteristic functional group is R–O–C(=O)–NR′R″.
"""
from rdkit import Chem

def is_carbamate_ester(smiles: str):
    """
    Determines if a molecule is a carbamate ester based on its SMILES string.
    A carbamate ester has the functional group R–O–C(=O)–NR', corresponding to an ester derived from carbamic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a carbamate ester, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for a carbamate ester:
    # [O;X2]   : An oxygen with two connections (ester oxygen)
    # C(=O)    : A carbonyl group
    # [NX3]    : A nitrogen with three connections (can be substituted)
    carbamate_pattern = Chem.MolFromSmarts("[O;X2]C(=O)[NX3]")
    if carbamate_pattern is None:
        return False, "Error in SMARTS pattern creation"
    
    # Check if the molecule contains the carbamate ester substructure
    if mol.HasSubstructMatch(carbamate_pattern):
        return True, "Contains carbamate ester functional group (R-O-C(=O)-NR')"
    else:
        return False, "Carbamate ester pattern not found in the molecule"
        
# Example usage (for testing purposes; these lines can be commented out or removed):
if __name__ == "__main__":
    test_smiles = [
        "CCOC(N)=O",  # urethane
        "COc1ccccc1OCC(O)COC(N)=O"  # 2-hydroxy-3-(2-methoxyphenoxy)propyl carbamate
    ]
    for smi in test_smiles:
        result, reason = is_carbamate_ester(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")