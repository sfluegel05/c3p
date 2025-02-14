"""
Classifies: CHEBI:35757 monocarboxylic acid anion
"""
from rdkit import Chem

def is_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a monocarboxylic acid anion based on its SMILES string.
    A monocarboxylic acid anion must have exactly one carboxylate group and no compensating
    positive charge if not a zwitterion.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a monocarboxylic acid anion, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Carboxylate group pattern
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=O)[O-]")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    
    # Ensure there's exactly one carboxylate and assess context
    if len(carboxylate_matches) == 1:
        # Check for other anionic groups not part of carboxylates
        extra_anionic_pattern = Chem.MolFromSmarts("[O-;!$([CX3](=O)[O-])]")        
        if mol.HasSubstructMatch(extra_anionic_pattern):
            return False, "Contains additional non-carboxylate anionic group"
        
        # Exclude zwitterions with positive charges
        zwitterion_pattern = Chem.MolFromSmarts("[N+]")
        if mol.HasSubstructMatch(zwitterion_pattern):
            return False, "Contains zwitterionic positive charges"
        
        return True, "Contains exactly one carboxylate group"
    
    # Multiple carboxylate groups might indicate more complex classes
    if len(carboxylate_matches) > 1:
        return False, "Found more than one carboxylate group"
    
    return False, "No carboxylate group found"