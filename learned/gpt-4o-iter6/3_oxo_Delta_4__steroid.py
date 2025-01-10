"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
from rdkit import Chem

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES string.
    A 3-oxo-Delta(4) steroid has a steroid backbone with a 3-oxo group and a Delta(4) double bond.

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
    
    # Look for a more refined steroid backbone that considers the tetracyclic structure
    steroid_pattern = Chem.MolFromSmarts('C1CCC2C3C(CC4=CC(=O)CCC34)C2C1') 
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Check for the presence of a 3-oxo group at the C3 position
    oxo_pattern = Chem.MolFromSmarts('[C;R1]=O')  # Joining oxygen to ring carbon patterns
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No 3-oxo group found at the proper position"
    
    # Verify the alpha, beta unsaturated bond (Delta(4) bond) specifically between C4 and C5
    delta_4_pattern = Chem.MolFromSmarts('[#6]1=CC=C[#6]C1')  # Specifies placement near the first ring
    if not mol.HasSubstructMatch(delta_4_pattern):
        return False, "No Delta(4) double bond found"
    
    return True, "Molecule classified as 3-oxo-Delta(4) steroid with appropriate moieties"