"""
Classifies: CHEBI:3098 bile acid
"""
"""
Classifies: bile acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    Bile acids are hydroxy-5beta-cholanic acids occurring in bile.
    This function checks for the presence of the steroid nucleus with
    correct stereochemistry (5β-configuration), hydroxyl groups, and
    a side chain ending with a carboxylic acid or its derivatives.

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

    # Define SMARTS pattern for steroid nucleus with 5β-configuration
    steroid_smarts = "[C@H]1CC[C@H]2[C@@H]([C@H]1)CC[C@H]3[C@@H]2CC[C@@]4(C)CC[C@@H](O)C4C3"
    steroid_pattern = Chem.MolFromSmarts(steroid_smarts)
    if steroid_pattern is None:
        return False, "Error parsing steroid SMARTS pattern"

    # Check for steroid nucleus with correct stereochemistry
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Steroid nucleus with 5β-configuration not found"

    # Check for carboxylic acid or derivatives in side chain
    # Carboxylic acid pattern
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[O;H]')
    # Ester pattern
    ester_pattern = Chem.MolFromSmarts('C(=O)O[C]')
    # Amide pattern (including glycine and taurine conjugates)
    amide_pattern = Chem.MolFromSmarts('C(=O)N')
    # Sulfonate pattern (for taurine conjugates)
    sulfonate_pattern = Chem.MolFromSmarts('S(=O)(=O)[O][C,N]')

    side_chain_match = False
    for pattern in [carboxylic_acid_pattern, ester_pattern, amide_pattern, sulfonate_pattern]:
        if mol.HasSubstructMatch(pattern):
            side_chain_match = True
            break

    if not side_chain_match:
        return False, "No carboxylic acid or derivative found in side chain"

    # Check for at least one hydroxyl group (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts('[OX2H]')
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    total_hydroxyls = len(hydroxyl_matches)
    if total_hydroxyls == 0:
        return False, "No hydroxyl groups found"

    return True, "Molecule matches bile acid structure with correct stereochemistry and side chain"

__metadata__ = {
    'chemical_class': {
        'name': 'bile acid',
        'definition': "Any member of a group of hydroxy-5beta-cholanic acids occurring in bile, where they are present as the sodium salts of their amides with glycine or taurine. In mammals bile acids almost invariably have 5beta-configuration."
    }
}