"""
Classifies: CHEBI:32863 secondary amine
"""
from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule is a secondary amine based on its SMILES string.
    A secondary amine is defined as having a nitrogen atom bonded to exactly two carbon atoms and one hydrogen atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define secondary amine pattern: [NH] connected to two carbon atoms
    secondary_amine_pattern = Chem.MolFromSmarts("[N;H1;R0]-[C]-[C]")
    
    if mol.HasSubstructMatch(secondary_amine_pattern):
        return True, "Contains a secondary amine group"
    else:
        return False, "Does not contain a secondary amine group"

# You can test the function with a SMILES string to see if it correctly classifies secondary amines.