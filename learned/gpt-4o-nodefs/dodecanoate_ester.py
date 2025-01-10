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
    
    # Look for an ester group pattern specifically connected to 12-carbons
    # Improved pattern: Dodecyl group connected to an ester (O=C-O-C12H25)
    dodecanoate_ester_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCC(=O)O")
    if not mol.HasSubstructMatch(dodecanoate_ester_pattern):
        return False, "No dodecanoic acid ester moiety found"
    
    return True, "Contains a dodecanoate ester moiety"

# Example usage
smiles_example = "CCCCCCCCCCCC(=O)OCC"
result, reason = is_dodecanoate_ester(smiles_example)
print(result, reason)