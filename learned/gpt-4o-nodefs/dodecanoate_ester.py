"""
Classifies: CHEBI:87659 dodecanoate ester
"""
from rdkit import Chem

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule is a dodecanoate ester based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dodecanoate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for an ester group pattern: -C(=O)O-
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester group found"
    
    # Check for a linear hydrocarbon chain of at least 12 carbons connected to the ester
    dodecanoic_acid_pattern = Chem.MolFromSmarts("CCCCCCCCCCCC[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(dodecanoic_acid_pattern):
        return False, "No dodecanoic acid moiety found"
    
    return True, "Contains an ester group with a dodecanoic acid moiety"

# Example usage
smiles_example = "CCCCCCCCCCCC(=O)OCC"
result, reason = is_dodecanoate_ester(smiles_example)
print(result, reason)