"""
Classifies: CHEBI:47788 3-oxo steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    A 3-oxo steroid has an oxo group at the 3-position on the steroid structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo steroid, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a general pattern for a steroid backbone (four fused rings with some flexibility)
    steroid_pattern = Chem.MolFromSmarts("C1CC2CCC3C4CCC(C=O)C4CCC3C2C=C1") # typical cyclopentanoperhydrophenanthrene structure
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Define a more flexible pattern for an oxo group at position 3 to cover diverse structures
    oxo_3_pattern = Chem.MolFromSmarts("C1=CC(=O)CC1") # Testing if 3-position has a carbonyl group connected in a ring system
    matches = mol.GetSubstructMatches(oxo_3_pattern)
    if not any(mol.GetAtomWithIdx(match[1]).GetNeighbors() and len(mol.GetAtomWithIdx(match[1]).GetNeighbors()) == 3 for match in matches):
        return False, "No 3-oxo group found on steroid skeleton"

    return True, "3-oxo group found at position 3 on steroid skeleton"

# End of code block