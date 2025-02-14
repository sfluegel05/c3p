"""
Classifies: CHEBI:35757 monocarboxylic acid anion
"""
from rdkit import Chem

def is_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a monocarboxylic acid anion based on its SMILES string.
    A monocarboxylic acid anion has exactly one carboxylate group, represented by the pattern C(=O)[O-],
    and lacks additional non-carboxylate charged functional groups.

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
    
    # Identify carboxylate group pattern with refinement for specific ion environments
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=O)[O-]")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    
    # Check for exactly one carboxylate group
    if len(carboxylate_matches) != 1:
        if len(carboxylate_matches) == 0:
            return False, "No carboxylate group found"
        else:
            return False, f"Found {len(carboxylate_matches)} carboxylate groups, expected exactly 1"
    
    # Identify and exclude additional anionic groups that aren't part of a carboxylate
    extra_anionic_pattern = Chem.MolFromSmarts("[O-;!$([CX3](=O)[O-])]")
    if mol.HasSubstructMatch(extra_anionic_pattern):
        return False, "Contains additional non-carboxylate anionic group indicating it's not a simple monocarboxylic acid anion"
    
    # Check for zwitterionic state by scanning common positive charge
    zwitterion_pattern = Chem.MolFromSmarts("[N+]")
    if mol.HasSubstructMatch(zwitterion_pattern):
        return False, "Contains zwitterionic structure indicating compensating positive charges exist"

    return True, "Contains exactly one carboxylate group, indicating a monocarboxylic acid anion"