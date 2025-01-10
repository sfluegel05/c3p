"""
Classifies: CHEBI:28863 flavanones
"""
from rdkit import Chem

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

    # Generalized chroman-4-one core pattern
    chromanone_pattern = Chem.MolFromSmarts("O=C1CCC2C1C=CC=C2") 
    if not mol.HasSubstructMatch(chromanone_pattern):
        return False, "No chroman-4-one core detected"
    
    # Phenyl group attached to the 2-position of the chroman-4-one
    phenyl_at_position_2_pattern = Chem.MolFromSmarts("c1ccccc1[C@H]2CC(O)=C2")
    if not mol.HasSubstructMatch(phenyl_at_position_2_pattern):
        return False, "No phenyl group at position 2 of the chroman-4-one core"

    return True, "Contains chroman-4-one core with a phenyl group at position 2"