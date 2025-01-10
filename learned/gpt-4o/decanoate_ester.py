"""
Classifies: CHEBI:87658 decanoate ester
"""
from rdkit import Chem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    A decanoate ester results from the esterification of decanoic acid with an alcohol or phenol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a decanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the decanoic acid moiety with variability in chain and possible connections
    decanoic_acid_patterns = [
        Chem.MolFromSmarts("CCCCCCCCCC(=O)O"),  # Simple linear
        Chem.MolFromSmarts("CCCCCCCC(CC)(=O)O"),  # Potential branching
        Chem.MolFromSmarts("C(=O)OCCCCCCCCCC")  # Reverse direction
    ]
    
    # Check if any of the decanoic acid patterns match
    decanoic_acid_found = False
    for pattern in decanoic_acid_patterns:
        if mol.HasSubstructMatch(pattern):
            decanoic_acid_found = True
            break
    
    if not decanoic_acid_found:
        return False, "No decanoic acid moiety found"

    # Ensure the ester linkage is linked to the decanoic acid moiety
    ester_linkage_pattern = Chem.MolFromSmarts("C(=O)O")  # This pattern might still be okay assuming it's matched in context

    if not mol.HasSubstructMatch(ester_linkage_pattern):
        return False, "No ester linkage found"

    return True, "Contains decanoic acid moiety with ester linkage"