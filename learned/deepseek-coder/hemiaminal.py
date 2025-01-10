"""
Classifies: CHEBI:73080 hemiaminal
"""
"""
Classifies: CHEBI:52217 hemiaminal
"""
from rdkit import Chem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule is a hemiaminal based on its SMILES string.
    A hemiaminal is any organic amino compound that has an amino group and a hydroxy group attached to the same carbon atom.

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

    # Define the hemiaminal pattern: a carbon atom bonded to exactly one amino group and one hydroxy group
    hemiaminal_pattern = Chem.MolFromSmarts("[C;H0;X4]([N;H1,H2])([O;H1])")
    
    # Check if the molecule contains the hemiaminal pattern
    if mol.HasSubstructMatch(hemiaminal_pattern):
        return True, "Contains a carbon atom bonded to both an amino group and a hydroxy group"
    else:
        return False, "No carbon atom found with both an amino group and a hydroxy group attached"