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

    # Improve pattern: Look for ester connected to a 12-carbon chain, taking into account variations
    dodecanoate_ester_pattern = Chem.MolFromSmarts("C(=O)O.CCCCCCCCCCCC")
    
    # Check for a match
    has_dodecanoate_ester = mol.HasSubstructMatch(dodecanoate_ester_pattern)
    if not has_dodecanoate_ester:
        return False, "No dodecanoic acid ester moiety found"
    
    return True, "Contains a dodecanoate ester moiety"

# Example usage
smiles_example = "CCCCCCCCCCCC(=O)OCC"
result, reason = is_dodecanoate_ester(smiles_example)
print(result, reason)