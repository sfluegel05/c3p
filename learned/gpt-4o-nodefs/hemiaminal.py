"""
Classifies: CHEBI:73080 hemiaminal
"""
from rdkit import Chem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule is a hemiaminal based on its SMILES string.
    A hemiaminal contains both a hydroxyl group (-OH) and an amine group (-NH2)
    attached to the same carbon.
    
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
    
    # Define hemiaminal pattern: C with OH and NH2 attached
    hemiaminal_pattern = Chem.MolFromSmarts("[C;R0]([OH])([NH2])")

    # Check for hemiaminal substructure
    if mol.HasSubstructMatch(hemiaminal_pattern):
        return True, "Contains hemiaminal substructure (-C(OH)(NH2)-)"
    else:
        return False, "Does not contain hemiaminal substructure"