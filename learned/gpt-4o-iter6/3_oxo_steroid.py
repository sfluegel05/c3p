"""
Classifies: CHEBI:47788 3-oxo steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

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
    
    # Define a more general pattern for the steroid backbone: four fused rings, allowing for variation
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C3C(C1)CCC3CCC4=C2C=CC4")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Define a pattern for an oxo group at position 3
    # This should involve a carbonyl bonded at the 3 position, recognizing typical labeling
    oxo_pattern = Chem.MolFromSmarts("C[C@@H]1CCCC2=C1C=CC(=O)C2")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No 3-oxo group found on steroid skeleton"

    return True, "3-oxo group found at position 3 on steroid skeleton"