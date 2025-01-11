"""
Classifies: CHEBI:51689 enone
"""
from rdkit import Chem

def is_enone(smiles: str):
    """
    Determines if a molecule is an enone based on its SMILES string.
    An enone is characterized by the presence of a C=O group adjacent to a C=C group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains an enone group, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the enone pattern: C=O adjacent to C=C
    enone_pattern = Chem.MolFromSmarts("C=CC(=O)")
    if mol.HasSubstructMatch(enone_pattern):
        return True, "Contains enone pattern: C=O adjacent to C=C"

    return False, "No enone pattern found"

# Example usage
example_smiles = "CC(C)[C@@H]1CC[C@H](C)[C@@]11CC=C(C)C(=O)C1"
is_enone(example_smiles)