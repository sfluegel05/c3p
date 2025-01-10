"""
Classifies: CHEBI:32863 secondary amine
"""
from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule is a secondary amine based on its SMILES string.
    A secondary amine is characterized by a nitrogen atom bonded to two carbon atoms and one hydrogen atom.

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
        
    # Look for secondary amine pattern: [NH](C)C
    secondary_amine_pattern = Chem.MolFromSmarts("[NH](C)C")
    if mol.HasSubstructMatch(secondary_amine_pattern):
        return True, "Contains a nitrogen atom bonded to two carbon atoms and one hydrogen atom, characteristic of a secondary amine"
    else:
        return False, "Does not satisfy the structural requirements for a secondary amine"