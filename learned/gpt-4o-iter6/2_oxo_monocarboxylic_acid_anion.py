"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid anion based on its SMILES string.
    This involves identifying a 2-oxo group at the correct position and a carboxylate anion.

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
    
    # Look for 2-oxo group adjacent to carboxylate (-C-C(=O)-C(=O)[O-])
    # Here we aim to have a carbon adjacent to the keto group 
    two_oxo_pattern = Chem.MolFromSmarts("[#6]-[#6](=O)-[#6](=O)[O-]")  # Specific C-C(=O)C(=O)[O-] arrangement
    if not mol.HasSubstructMatch(two_oxo_pattern):
        return False, "2-oxo group not found at expected position relative to carboxylate"

    return True, "Contains a 2-oxo group at the 2-position and a carboxylate anion"

# Examples for testing
# should return True for molecules with 2-oxo monocarboxylic acid anions
#print(is_2_oxo_monocarboxylic_acid_anion('CNC(=O)CCC(=O)C([O-])=O')) # True
#print(is_2_oxo_monocarboxylic_acid_anion('CC(=O)C(C)(O)C([O-])=O')) # True
#print(is_2_oxo_monocarboxylic_acid_anion('CC(C)CC(=O)O')) # False