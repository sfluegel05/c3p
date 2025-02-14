"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.
    A 2-oxo monocarboxylic acid contains one carboxylic acid group and a carbonyl group at
    the alpha position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-oxo monocarboxylic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for monocarboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) != 1:
         return False, f"Molecule must have exactly one carboxylic acid group, found {len(carboxylic_acid_matches)}"

    # 2. Check for 2-oxo substituent (carbonyl at alpha position) - using smarts to find chain
    oxo_group_pattern = Chem.MolFromSmarts("C(=O)C(=O)O")
    oxo_matches = mol.GetSubstructMatches(oxo_group_pattern)
    if len(oxo_matches) < 1:
        # check for carbonyl immediately adjacent to carboxyl group
        oxo_pattern_2 = Chem.MolFromSmarts("[CX3](=O)[CX4][CX3](=O)[OX2]")
        oxo_matches_2 = mol.GetSubstructMatches(oxo_pattern_2)
        if len(oxo_matches_2) == 0:
            return False, "No 2-oxo group found"


    return True, "Molecule has a 2-oxo substituent and one carboxylic acid group"