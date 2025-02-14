"""
Classifies: CHEBI:26004 phenylpropanoid
"""
from rdkit import Chem

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    A phenylpropanoid typically has a phenylpropane backbone with various functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a phenylpropanoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for a benzene ring
    aromatic_ring_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(aromatic_ring_pattern):
        return False, "No aromatic benzene ring found"

    # Look for a propane chain attached to the aromatic ring
    phenylpropane_pattern_1 = Chem.MolFromSmarts("c1ccccc1CCC")
    phenylpropane_pattern_2 = Chem.MolFromSmarts("c1ccccc1[CH2,CH][CH2,CH]")  # More flexible version
    if not (mol.HasSubstructMatch(phenylpropane_pattern_1) or mol.HasSubstructMatch(phenylpropane_pattern_2)):
        return False, "No phenylpropane backbone found"

    # Check for hydroxyl, methoxy or other functional groups commonly found in phenylpropanoids
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    methoxy_pattern = Chem.MolFromSmarts("CO")
    
    contains_functional_groups = (
        mol.HasSubstructMatch(hydroxyl_pattern) or mol.HasSubstructMatch(methoxy_pattern)
    )
    
    if not contains_functional_groups:
        return False, "Missing characteristic functional groups (hydroxyl/methoxy)"

    return True, "Contains a phenylpropane backbone and characteristic functional groups"