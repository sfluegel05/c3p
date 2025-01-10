"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid anion based on its SMILES string.
    This involves identifying a 2-oxo group and a carboxylate anion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-oxo monocarboxylic acid anion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylate anion (-C(=O)[O-])
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate anion found"
        
    # Look for 2-oxo ketone (C-C(=O)-) next to carboxylate
    # Constructing a pattern that seeks a sequence of 'any atom, C, and O double bonded' 
    # adjacent to another carbon, representing '-C-C(=O)-C(=O)'
    two_oxo_pattern = Chem.MolFromSmarts("C-C(=O)C(=O)[O-]")
    if not mol.HasSubstructMatch(two_oxo_pattern):
        return False, "2-oxo group not found at expected position"
    
    return True, "Contains a 2-oxo group and a carboxylate anion"

# Examples for testing
# print(is_2_oxo_monocarboxylic_acid_anion('CNC(=O)CCC(=O)C([O-])=O')) # Should return True
# print(is_2_oxo_monocarboxylic_acid_anion('CC(C)CC(=O)O')) # Should return False, not an anion