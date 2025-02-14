"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
from rdkit import Chem

def is_3_oxo_5alpha_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5alpha-steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a 3-oxo-5alpha-steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string into molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # General steroid backbone: four-ring structure (A, B, C, D)
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C3C4C(C)CCC4C3CCC2C1")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Check for 3-oxo group ==O on the C attached to second ring (consistent with steroid A ring)
    # This pattern looks for a carbonyl group on C3
    oxo_pattern = Chem.MolFromSmarts("C(=O)[C@H]")
    oxo_matches = mol.GetSubstructMatch(oxo_pattern)
    if not oxo_matches or len(oxo_matches) != 1 or mol.GetAtomWithIdx(oxo_matches[0]).GetIsAromatic():
        return False, "No 3-oxo group found"

    # Validate 5alpha configuration
    # The 5alpha configuration requires a hydrogen to be in the 'alpha' position
    five_alpha_pattern = Chem.MolFromSmarts("[C@@H](C)C1CCC2")
    if not mol.HasSubstructMatch(five_alpha_pattern):
        return False, "No 5alpha configuration found"

    return True, "Contains a 3-oxo group and the 5alpha stereochemical configuration, typical of 3-oxo-5alpha-steroids"