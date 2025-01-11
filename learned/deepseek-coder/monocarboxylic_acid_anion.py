"""
Classifies: CHEBI:35757 monocarboxylic acid anion
"""
"""
Classifies: CHEBI:51086 monocarboxylic acid anion
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

    # Check for the presence of a carboxylate group (C(=O)[O-])
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[O-]")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    
    if len(carboxylate_matches) != 1:
        return False, f"Found {len(carboxylate_matches)} carboxylate groups, need exactly 1"

    # Ensure there are no other negatively charged atoms
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() == -1 and not atom.GetSymbol() == 'O':
            return False, "Found negatively charged atom that is not part of the carboxylate group"

    # Ensure there are no other carboxyl groups (COOH)
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    
    if len(carboxyl_matches) > 0:
        return False, "Found additional carboxyl groups (COOH)"

    return True, "Contains exactly one carboxylate group and no other negatively charged atoms"