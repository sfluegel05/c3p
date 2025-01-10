"""
Classifies: CHEBI:73080 hemiaminal
"""
from rdkit import Chem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule is a hemiaminal based on its SMILES string.
    A hemiaminal is characterized by an amino group and a hydroxy group attached to the same carbon atom.

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

    # Define hemiaminal pattern: Carbon attached to both nitrogen (N) and oxygen (O)
    hemiaminal_pattern = Chem.MolFromSmarts("[CX4]([NH2,NH,N])([OH])")
    if not mol.HasSubstructMatch(hemiaminal_pattern):
        return False, "No hemiaminal structure found: missing amino and hydroxy groups on the same carbon"

    return True, "Contains hemiaminal structure: amino and hydroxy groups attached to the same carbon"