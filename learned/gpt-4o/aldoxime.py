"""
Classifies: CHEBI:22307 aldoxime
"""
from rdkit import Chem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime based on its SMILES string.
    An aldoxime is an oxime of an aldehyde containing the group RCH=NOH.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldoxime, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Updated pattern for an aldoxime focuses on the presence of RCH=NOH
    aldoxime_pattern = Chem.MolFromSmarts("[#6][CX3](=N[OH1])[#1]") 

    if mol.HasSubstructMatch(aldoxime_pattern):
        return True, "Contains C=NOH group indicating aldoxime structure"

    return False, "Does not contain the C=NOH group typical of an aldoxime"

# Sample test
smiles = "C(C(C)C)=NO"  # example of an aldoxime
result, reason = is_aldoxime(smiles)
print(result, reason)