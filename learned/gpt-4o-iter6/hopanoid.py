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

    # Enhanced hopanoid patterns capturing more comprehensive diversity
    hopanoid_patterns = [
        Chem.MolFromSmarts('[C@@]12CC[C@H]3[C@@]([C@H]4CC[C@@]5([C@@H](CC[C@]45)C)C[C@@]3(C)CC2)C1'),
        Chem.MolFromSmarts('[C@]12CCC[C@@]3([C@@H](CC4([C@@H]([C@H](CC34)C)C)C[C@@]2C)C1)C'),
        # Include patterns with common functional groups or modifications
        Chem.MolFromSmarts('[C@@]12CC[C@H]3[C@]4([C@@H](C[C@H](O)[C@@]45)C)CCC[C@]35(C)CC2)C1'),
    ]

    # Check if any of the patterns match
    for pattern in hopanoid_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains hopanoid hopane skeleton based on pentacyclic triterpenoid structure"

    return False, "Hopanoid hopane skeleton not found"