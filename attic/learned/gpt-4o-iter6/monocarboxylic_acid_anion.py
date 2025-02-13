"""
Classifies: CHEBI:35757 monocarboxylic acid anion
"""
from rdkit import Chem

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
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for carboxylate group [O-]C(=O)
    carboxylate_pattern = Chem.MolFromSmarts("[O-]C(=O)")
    
    # Find matches for carboxylate pattern
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    if len(carboxylate_matches) == 0:
        return False, "No carboxylate group found"
    if len(carboxylate_matches) > 1:
        return False, f"Multiple ({len(carboxylate_matches)}) carboxylate groups found"

    return True, "Contains a single deprotonated carboxylic acid group"