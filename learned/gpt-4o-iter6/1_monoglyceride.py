"""
Classifies: CHEBI:35759 1-monoglyceride
"""
from rdkit import Chem

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 1-monoglyceride based on its SMILES string.
    A 1-monoglyceride is defined as a monoglyceride in which the acyl substituent is located at position 1.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-monoglyceride, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for extracting the 1-monoglyceride features:
    # - Primary carbon ester linkage [C@](O)COC(=O)
    # - Secondary alcohol remains [C@@](O)COH
    # - Exact chain length or variations

    # Look for 1-monoglyceride specific structural patterns
    monoglyceride_pattern = Chem.MolFromSmarts('O[C@@H](CO)COC(=O)C')  # Active ester and secondary alcohol pattern
    if not mol.HasSubstructMatch(monoglyceride_pattern):
        return False, "No 1-monoglyceride backbone detected"

    # Validity check on the specific esterified position:
    ester_pattern = Chem.MolFromSmarts('OC([C@H](CO)C(=O)O)')
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "Acyl linkage not located at position 1"

    # Check for additional branches if any so that primary criteria are met
    branches = mol.GetSubstructMatches(ester_pattern)
    if len(branches) != 1:
        return False, f"Found {len(branches)} branches of esterification, expected exactly 1 at the primary position"

    return True, "Verified 1-monoglyceride with acyl linkage at position 1"