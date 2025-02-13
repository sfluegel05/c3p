"""
Classifies: CHEBI:22307 aldoxime
"""
from rdkit import Chem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime based on its SMILES string.
    An aldoxime is an oxime of an aldehyde RCH=NOH.

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

    # The optimized pattern reflects an aldoxime functional group focusing on C=N-O-H 
    # where C must be attached directly to only one other atom (single-bonded) or hydrogen to reflect CH=NOH
    aldoxime_pattern = Chem.MolFromSmarts("[#6][CX2,H]=[NX2]([OX1H])") 

    if mol.HasSubstructMatch(aldoxime_pattern):
        return True, "Contains C=NOH group indicating aldoxime structure"

    return False, "Does not contain the C=NOH group typical of an aldoxime"

# Sample test
smiles = "C(C(C)C)=NO"  # example of an aldoxime
result, reason = is_aldoxime(smiles)
print(result, reason)