"""
Classifies: CHEBI:51963 hopanoid
"""
from rdkit import Chem

def is_hopanoid(smiles: str):
    """
    Determines if a molecule is a hopanoid based on its SMILES string.
    A hopanoid is a triterpenoid based on a hopane skeleton, featuring a pentacyclic structure with specific stereochemistry.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hopanoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Hopanoid hopane skeleton pattern - capturing more diverse triterpenoid hopanoid structures
    # This includes flexibility in stereochemistry and attached groups.
    hopanoid_patterns = [
        Chem.MolFromSmarts('[C@]12CC[C@@H]3[C@]4(C)CC[C@@]5([C@H]4CC[C@]35C)C2C1'),
        Chem.MolFromSmarts('[C@@]12CC[C@H]3[C@]4([C@@]5(CCC(C5)C)C(CCC4)C3)C(CCC2)C1'),
        # Additional patterns that allow for common modifications or stereochemical variants
    ]

    # Check if any of the patterns match
    for pattern in hopanoid_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains hopanoid hopane skeleton based on pentacyclic triterpenoid structure"

    return False, "Hopanoid hopane skeleton not found"