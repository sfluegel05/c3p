"""
Classifies: CHEBI:23003 carbamate ester
"""
from rdkit import Chem

def is_carbamate_ester(smiles: str):
    """
    Determines if a molecule is a carbamate ester based on its SMILES string.
    A carbamate ester contains the functional group -OC(=O)N- which is an ester 
    of carbamic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a carbamate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a carbamate ester
    carbamate_pattern = Chem.MolFromSmarts('OC(=O)N')
    
    # Check if the molecule contains the carbamate pattern
    if mol.HasSubstructMatch(carbamate_pattern):
        return True, "Molecule contains the carbamate ester functional group"
    else:
        return False, "No carbamate ester functional group found"

# Example test cases - Uncomment to test
# print(is_carbamate_ester('CNC(=O)Oc1ccc(C)c(C)c1')) # Should return True
# print(is_carbamate_ester('CC(=O)OC1=CC=CC=C1C(=O)O')) # Should return False