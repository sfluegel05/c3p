"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
from rdkit import Chem

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid anion based on its SMILES string.
    This involves identifying a 2-oxo group at the 2-position and a carboxylate anion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 2-oxo monocarboxylic acid anion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylate anion (-C(=O)[O-]) pattern
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate anion found"
    
    # Look for 2-oxo group pattern adjacent to carboxylate (-[CH2,C](=O)C(=O)[O-])
    # Allow flexibility in the carbon backbone, detecting any potential substitutions around the C=O group
    two_oxo_pattern = Chem.MolFromSmarts("[#6]-[#6X3v4H0,H1,H2]-C(=O)-[#6](=O)[O-]")  
    if not mol.HasSubstructMatch(two_oxo_pattern):
        return False, "2-oxo group not at expected position relative to carboxylate"

    return True, "Contains a 2-oxo group at the 2-position and a carboxylate anion"

# Test cases for verification can be added by the user in a testing script or environment