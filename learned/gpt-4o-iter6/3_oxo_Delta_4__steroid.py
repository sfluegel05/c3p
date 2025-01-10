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
    
    # Define a specific steroid backbone pattern for the tetracyclic ring system
    steroid_pattern = Chem.MolFromSmarts('C1CCC2C3CCC4CCCC(C3)C2C1C4')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Check for the presence of a 3-oxo group specifically on the A-ring
    oxo_pattern = Chem.MolFromSmarts('C1(=O)[C,C]C2')
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No 3-oxo group found in the appropriate position"
    
    # Verify the Alpha, Beta unsaturated bond (Delta(4) bond)
    # This means C=C between carbons 4 and 5 in the steroid A-ring
    delta_4_pattern = Chem.MolFromSmarts('C1=C[C,C]2')
    if not mol.HasSubstructMatch(delta_4_pattern):
        return False, "No Delta(4) double bond found in the appropriate position"
    
    return True, "Molecule classified as 3-oxo-Delta(4) steroid with appropriate moieties"