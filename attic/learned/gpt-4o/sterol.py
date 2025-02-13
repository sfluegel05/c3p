"""
Classifies: CHEBI:15889 sterol
"""
from rdkit import Chem

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.
    A sterol is a 3-hydroxy steroid with a structure related to cholestan-3-ol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for steroid backbone (the rings ABCD)
    steroid_pattern = Chem.MolFromSmarts('[#6]1[#6][#6][#6]2[#6][#6]3[#6][#6][#6]4[#6][#6][#6][#6]4[#6]3[#6][#6][#6]21')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # SMARTS pattern for hydroxyl group on the A-ring
    hydroxyl_pattern = Chem.MolFromSmarts('[#6][#6]([O])')  # Assumes hydroxyl group is attached to a carbon
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No 3-hydroxy group found"

    # Check for potential side chains (not strict here, but checking continuity of side chain)
    side_chain = any(atom.GetDegree() > 3 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not side_chain:
        return False, "No suitable side chain found"

    return True, "Structure closely related to a sterol"

# Examples given can be tested with this function to check the classification