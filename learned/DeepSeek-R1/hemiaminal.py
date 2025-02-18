"""
Classifies: CHEBI:73080 hemiaminal
"""
"""
Classifies: CHEBI:50047 hemiaminal
"""
from rdkit import Chem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule is a hemiaminal based on its SMILES string.
    A hemiaminal has an amino group and a hydroxy group attached to the same carbon atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hemiaminal, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for a carbon with both hydroxyl and amino groups
    # [C] connected to [OH] and [NX3] (amines: NH2, NHR, NR2)
    hemiaminal_pattern = Chem.MolFromSmarts("[C]([OH])([NX3])")
    
    if mol.HasSubstructMatch(hemiaminal_pattern):
        return True, "Contains a carbon with both hydroxyl and amino groups"
    else:
        return False, "No carbon with hydroxyl and amino groups found"