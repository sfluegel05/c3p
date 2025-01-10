"""
Classifies: CHEBI:73080 hemiaminal
"""
from rdkit import Chem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule is a hemiaminal based on its SMILES string.

    A hemiaminal contains both a hydroxyl group (-OH) and an amine group (-NH2 or -NH)
    attached to the same carbon atom, considering possible stereochemistry.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hemiaminal, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Unified pattern for hemiaminal:
    # Generalized pattern accommodating stereochemistry and diverse structural contexts
    hemiaminal_pattern = Chem.MolFromSmarts("[C]([OH])([NH2,NH])")

    # Attempt to match the hemiaminal pattern
    if mol.HasSubstructMatch(hemiaminal_pattern):
        return True, "Contains detectable hemiaminal substructure"

    return False, "Does not contain detectable hemiaminal substructure"