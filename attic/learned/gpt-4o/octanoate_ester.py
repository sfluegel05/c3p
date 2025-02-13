"""
Classifies: CHEBI:87657 octanoate ester
"""
from rdkit import Chem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester contains the ester linkage where the carboxylic acid component is octanoic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define octanoic acid ester pattern: R-O-C(=O)-CCCCCCC
    octanoate_ester_pattern = Chem.MolFromSmarts("C(=O)OCCCCCCC")
    if mol.HasSubstructMatch(octanoate_ester_pattern):
        return True, "Contains octanoic acid ester linkage"
    
    return False, "Does not contain octanoic acid ester linkage"

# Test the function with sample SMILES strings
test_smiles = [
    "CCCCCCCC(=O)OC(COC(=O)CCCCCCC)OC(=O)CCCCCCC",  # Trioctanoin
    "CCCCCCCC(=O)OC",  # Methyl octanoate
    "CCCCCCCC",  # Octane, not an ester
]

for smiles in test_smiles:
    result, reason = is_octanoate_ester(smiles)
    print(f"SMILES: {smiles} -> Is Octanoate Ester: {result}, Reason: {reason}")