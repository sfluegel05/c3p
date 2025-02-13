"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
from rdkit import Chem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    A tetrahydrofuranone is any oxolane with an oxo-substituent on the tetrahydrofuran ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrahydrofuranone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Define flexible patterns for tetrahydrofuranone
    # Generic five-membered ring with one oxygen and one oxo group (anywhere on the ring)
    thf_patterns = [
        Chem.MolFromSmarts("C1C(=O)COC1"),    # Oxo on 2nd position
        Chem.MolFromSmarts("C1CC(=O)OC1"),    # Oxo on 3rd position
        Chem.MolFromSmarts("C1CCOC(=O)1"),    # Oxo on 4th position
        Chem.MolFromSmarts("C1(C=O)CCOC1"),   # Oxo on 5th position
    ]

    # Ensure patterns are valid
    if not all(p for p in thf_patterns):
        return None, "Error in generating SMARTS pattern"

    # Check against each pattern
    for thf_pattern in thf_patterns:
        if mol.HasSubstructMatch(thf_pattern):
            return True, "Contains a tetrahydrofuran ring with an oxo group"

    return False, "Does not contain a tetrahydrofuran ring with an oxo group"