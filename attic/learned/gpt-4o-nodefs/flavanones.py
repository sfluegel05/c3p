"""
Classifies: CHEBI:28863 flavanones
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_flavanones(smiles: str):
    """
    Determines if a molecule is a flavanone based on its SMILES string.
    A flavanone has a chroman-4-one core with a phenyl group at position 2.

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

    # Pattern for the chroman-4-one core
    chromanone_pattern = Chem.MolFromSmarts("O=C1CC(OC2=CC=CC=C2)C=3C=CC(O)=CC3C1")
    if not mol.HasSubstructMatch(chromanone_pattern):
        return False, "No chroman-4-one core detected"
    
    # Check for the presence of a carbonyl group at position 4
    carbonyl_pattern = Chem.MolFromSmarts("O=CC2C=3C=CC(O)=CC3C2")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "Missing ketone at position 4"
    
    # Match a phenyl group at position 2
    phenyl_group_pattern = Chem.MolFromSmarts("C1=CC=CC=C1")
    phenyl_matches = mol.GetSubstructMatches(phenyl_group_pattern)
    if len(phenyl_matches) < 1:
        return False, "No phenyl group at position 2"

    return True, "Contains chroman-4-one core with a phenyl group at position 2"