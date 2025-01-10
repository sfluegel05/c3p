"""
Classifies: CHEBI:32876 tertiary amine
"""
from rdkit import Chem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is a tertiary amine based on its SMILES string.
    A tertiary amine has a nitrogen atom bonded to three carbon atoms,
    potentially including different types of carbon hybridizations.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Flexible pattern for tertiary amine: nitrogen with three carbon attachments,
    # allowing aliphatic, aromatic, and various hybridizations 
    tertiary_amine_pattern = Chem.MolFromSmarts("[NX3]([C,N,O,S,P])([C,N,O,S,P])([C,N,O,S,P])")
    if mol.HasSubstructMatch(tertiary_amine_pattern):
        return True, "Contains a nitrogen atom bonded to three groups, indicating a tertiary amine"

    # If no match is found
    return False, "Does not contain a tertiary amine group"