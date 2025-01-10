"""
Classifies: CHEBI:3098 bile acid
"""
"""
Classifies: bile acid
"""
from rdkit import Chem

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    Bile acids are hydroxy-5beta-cholanic acids occurring in bile.
    This function checks for the presence of a steroid nucleus with a carboxylic acid side chain and hydroxyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for steroid nucleus
    steroid_smarts = '[#6]1[#6][#6][#6][#6][#6]1[#6]2[#6][#6][#6][#6][#6]2[#6]3[#6][#6][#6][#6][#6]3[#6]4[#6][#6][#6][#6]4'
    steroid_pattern = Chem.MolFromSmarts(steroid_smarts)
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Steroid nucleus not found"

    # Check for carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[O;H]')
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Carboxylic acid group not found"

    # Check for hydroxyl groups (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts('[OX2H]')
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) == 0:
        return False, "No hydroxyl groups found"

    # Optionally, check for 5beta-configuration (complex stereochemistry)
    # Skipping detailed stereochemistry check for simplicity

    return True, "Molecule contains steroid nucleus with carboxylic acid and hydroxyl groups"


__metadata__ = {
    'chemical_class': {
        'name': 'bile acid',
        'definition': "Any member of a group of hydroxy-5beta-cholanic acids occurring in bile, where they are present as the sodium salts of their amides with glycine or taurine. In mammals bile acids almost invariably have 5beta-configuration."
    }
}