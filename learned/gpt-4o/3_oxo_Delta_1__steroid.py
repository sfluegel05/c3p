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

    # Improved SMARTS pattern for 3-oxo group (C=O on the steroid backbone)
    oxo_pattern = Chem.MolFromSmarts("C(=O)[C;R3]")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No 3-oxo group (ketone) found"

    # Improved SMARTS pattern for Delta(1) double bond
    double_bond_pattern = Chem.MolFromSmarts("[C;R1]=[C;R2]")
    if not mol.HasSubstructMatch(double_bond_pattern):
        return False, "No Delta(1) double bond (C=C between position 1 and 2) found"

    # Steroid core structure involves a fused four-ring system (tetracyclic)
    steroid_core_pattern = Chem.MolFromSmarts("C1CC2CCC3C(=O)CCC3C2C1")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Molecule does not have the typical steroid four-ring core structure"

    return True, "Matches 3-oxo-Delta(1) steroid structure"