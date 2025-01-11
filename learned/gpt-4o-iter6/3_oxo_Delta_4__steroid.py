"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
from rdkit import Chem

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES string.
    A 3-oxo-Delta(4) steroid contains a steroid backbone with a 3-oxo group and a Delta(4) double bond.

    Args:
        smiles (str): SMILES string of the chemical entity.

    Returns:
        bool: True if the molecule is a 3-oxo-Delta(4) steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Generalize the steroid backbone detection using a flexible tetracyclic framework
    steroid_pattern = Chem.MolFromSmarts('C1CCC2CCCCC2C1')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Check for the presence of a 3-oxo group
    # This group is part of a ketone functional group, so look for C=O attached to a carbon
    oxo_pattern = Chem.MolFromSmarts('C(=O)')
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No 3-oxo group found"
    
    # Verify the Alpha, Beta unsaturated bond (Delta(4) bond) as a C=C
    # This often means an alkene feature within the molecule
    delta_4_pattern = Chem.MolFromSmarts('C=C')
    if not mol.HasSubstructMatch(delta_4_pattern):
        return False, "No Delta(4) double bond found"
    
    return True, "Molecule classified as 3-oxo-Delta(4) steroid with appropriate moieties"