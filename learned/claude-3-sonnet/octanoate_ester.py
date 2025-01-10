"""
Classifies: CHEBI:87657 octanoate ester
"""
"""
Classifies: CHEBI:52247 octanoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester is any ester where the acid component is octanoic acid (caprylic acid).

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

    # Look for ester group pattern
    ester_pattern = Chem.MolFromSmarts("[#6][CX3](=[OX1])[OX2][#6]")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester group found"
    
    # Pattern for octanoate group: linear 7-carbon chain attached to ester carbonyl
    # Note: The SMARTS below looks for:
    # - An ester group (-C(=O)O-)
    # - Connected to exactly 6 more carbons in a linear chain
    # - The last carbon must have 3 hydrogens (CH3)
    octanoate_pattern = Chem.MolFromSmarts("[CH3][CH2][CH2][CH2][CH2][CH2][CH2]C(=O)[OX2]")
    
    if not mol.HasSubstructMatch(octanoate_pattern):
        return False, "No octanoate (8-carbon) chain found attached to ester group"
    
    # Count matches to ensure we have at least one octanoate ester
    matches = mol.GetSubstructMatches(octanoate_pattern)
    
    # Get the atoms involved in each match to verify the chain length
    for match in matches:
        # The pattern above ensures we have the right structure
        return True, "Contains octanoate ester group (8-carbon chain with terminal ester)"
        
    return False, "No valid octanoate ester found"

def test_examples():
    """Test function with some known examples"""
    test_cases = [
        ("CCCCCCCC(=O)OCC", True),  # ethyl octanoate
        ("CCCCCCCC(=O)OC", True),   # methyl octanoate
        ("CCCCCCC(=O)OCC", False),  # ethyl heptanoate (too short)
        ("CCCCCCCCC(=O)OCC", False),# ethyl nonanoate (too long)
        ("CCCCCCCCOC", False),      # octyl ether (not an ester)
    ]
    
    for smiles, expected in test_cases:
        result, reason = is_octanoate_ester(smiles)
        print(f"SMILES: {smiles}")
        print(f"Expected: {expected}, Got: {result}")
        print(f"Reason: {reason}\n")

if __name__ == "__main__":
    test_examples()