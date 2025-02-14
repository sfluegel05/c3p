"""
Classifies: CHEBI:20156 3-oxo-Delta(1) steroid
"""
from rdkit import Chem

def is_3_oxo_Delta_1__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(1) steroid based on its SMILES string.
    A 3-oxo-Delta(1) steroid is a steroid with a ketone group at position 3 and a double bond
    between positions 1 and 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-Delta(1) steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Updated SMARTS pattern for a 3-oxo group (carbonyl C=O at position 3)
    oxo_pattern = Chem.MolFromSmarts("C1=CC[C@]2(C)CC(=O)CC3=C2C[CH2]CC18")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No 3-oxo group (ketone) found on the third carbon"
    
    # Updated SMARTS pattern for a Delta(1) double bond (C=C between the first two carbons)
    double_bond_pattern = Chem.MolFromSmarts("C=C-C[CH2]")
    if not mol.HasSubstructMatch(double_bond_pattern):
        return False, "No Delta(1) double bond (C=C between C1 and C2) found"
    
    # Refined SMARTS pattern for the core tetracyclic structure of steroids
    steroid_core_pattern = Chem.MolFromSmarts("C1CCC2C3CCC4CCCC3C24C1")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Molecule does not have the steroid four-ring core structure"
    
    return True, "Matches 3-oxo-Delta(1) steroid structure"