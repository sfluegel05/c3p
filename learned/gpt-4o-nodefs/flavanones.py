"""
Classifies: CHEBI:28863 flavanones
"""
from rdkit import Chem

def is_flavanones(smiles: str):
    """
    Determines if a molecule is a flavanone based on its SMILES string.
    A flavanone features a 2,3-dihydro-2-phenylchromen-4-one core structure,
    often with various substitutions and stereochemistry.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavanone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for 2,3-dihydro-2-phenylchromen-4-one core
    chromanone_pattern = Chem.MolFromSmarts("O=C1C=C(OCC2=CC=CC=C2)C3=CC=CC=C3C1")

    # Check if molecule has the core chroman-4-one structure
    if not mol.HasSubstructMatch(chromanone_pattern):
        return False, "No recognized flavanone core detected"

    return True, "Contains recognizable flavanone core structure"