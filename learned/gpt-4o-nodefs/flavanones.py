"""
Classifies: CHEBI:28863 flavanones
"""
from rdkit import Chem

def is_flavanones(smiles: str):
    """
    Determines if a molecule is a flavanone based on its SMILES string.
    A flavanone typically features a chroman-4-one core structure,
    often with a phenyl or a similarly structured group.

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

    # Updated chroman-4-one core pattern to be more flexible
    chromanone_pattern = Chem.MolFromSmarts("O=C1CC[C@@]2(c3ccccc3)OCCc12") 
    if not mol.HasSubstructMatch(chromanone_pattern):
        return False, "No recognized chroman-4-one core detected"
    
    # Check for the phenyl-like ring attached to the core
    phenyl_pattern = Chem.MolFromSmarts("c1cc(ccc1)-C2OC(C)=C2")
    if not mol.HasSubstructMatch(phenyl_pattern):
        return False, "No phenyl or similar group detected at flavanone core"

    return True, "Contains recognizable flavanone core with appropriate substitution"