"""
Classifies: CHEBI:32876 tertiary amine
"""
from rdkit import Chem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is a tertiary amine based on its SMILES string.
    A tertiary amine is characterized by a nitrogen atom bonded specifically to three carbon atoms, 
    ignoring the context such as whether it is part of an amide or linked to sp2 hybrids.

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

    # Define a SMARTS pattern for tertiary amines: any nitrogen with three carbon bonds
    # Allow context of aromatic or amide connections as long as N has exactly three carbon attachments
    tertiary_amine_pattern = Chem.MolFromSmarts("[NX3]([CX3,CX4])([CX3,CX4])([CX3,CX4])")
    if mol.HasSubstructMatch(tertiary_amine_pattern):
        return True, "Contains a nitrogen atom bonded to three carbon atoms, indicating a tertiary amine"
    
    # If no match is found
    return False, "Does not contain a tertiary amine group"