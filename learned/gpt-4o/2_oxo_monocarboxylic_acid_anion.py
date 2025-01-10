"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
from rdkit import Chem

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid anion based on its SMILES string.
    This class has an oxo group located at the 2-position and a carboxylate group as the anion.

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

    # Look for a specific pattern that represents:
    # - a 2-oxo group pattern: carbonyl group bonded to another carbon at the 2-position
    # - followed by a variable hydrocarbon chain
    # - ending with a carboxylate group (correcting for adjacency is key)
    oxo_carboxylate_pattern = Chem.MolFromSmarts("C(=O)[CX3][CX4,C](=O)[O-]")
    
    if not mol.HasSubstructMatch(oxo_carboxylate_pattern):
        return False, "Pattern for 2-oxo monocarboxylic acid anion not matched"

    # Additional filtering criteria if necessary can be applied here
    
    return True, "Contains 2-oxo group and carboxylate anion in the correct context for a 2-oxo monocarboxylic acid anion"