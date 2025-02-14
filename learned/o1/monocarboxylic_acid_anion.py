"""
Classifies: CHEBI:35757 monocarboxylic acid anion
"""
"""
Classifies: monocarboxylic acid anion
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a monocarboxylic acid anion based on its SMILES string.
    A monocarboxylic acid anion is formed when the carboxy group of a monocarboxylic acid is deprotonated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monocarboxylic acid anion, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the pattern for any carboxyl group (protonated or deprotonated)
    carboxyl_group_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1H0-,OX1H1]")
    # Find all carboxyl groups in the molecule
    carboxyl_groups = mol.GetSubstructMatches(carboxyl_group_pattern)
    num_carboxyl_groups = len(carboxyl_groups)
    
    # Check that there is exactly one carboxyl group
    if num_carboxyl_groups != 1:
        return False, f"Found {num_carboxyl_groups} carboxyl groups, expected exactly 1"
    
    # Define the pattern for a carboxylate group (deprotonated carboxyl group)
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=O)[O-]")
    # Find all carboxylate groups in the molecule
    carboxylate_groups = mol.GetSubstructMatches(carboxylate_pattern)
    num_carboxylate_groups = len(carboxylate_groups)
    
    # Check that there is exactly one carboxylate group
    if num_carboxylate_groups != 1:
        return False, f"Found {num_carboxylate_groups} carboxylate groups, expected exactly 1"
    
    # If both conditions are met, the molecule is a monocarboxylic acid anion
    return True, "Molecule is a monocarboxylic acid anion with exactly one deprotonated carboxyl group"