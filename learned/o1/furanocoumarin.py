"""
Classifies: CHEBI:24128 furanocoumarin
"""
"""
Classifies: furanocoumarin
"""

from rdkit import Chem
from rdkit.Chem import rdqueries

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    A furanocoumarin consists of a furan ring fused with a coumarin. The fusion
    may occur in different ways, leading to several isomers (linear and angular).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a furanocoumarin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for linear and angular furanocoumarins
    # Linear furanocoumarin (psoralen-type) core SMARTS
    linear_smarts = 'c1cc2oc3ccoc3cc2c1'  # Matches psoralen core
    linear_pattern = Chem.MolFromSmarts(linear_smarts)

    # Angular furanocoumarin (angelicin-type) core SMARTS
    angular_smarts = 'c1ccc2oc3ccocc3c2c1'  # Matches angelicin core
    angular_pattern = Chem.MolFromSmarts(angular_smarts)

    if linear_pattern is None or angular_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Check for linear furanocoumarin core
    if mol.HasSubstructMatch(linear_pattern):
        return True, "Contains linear furanocoumarin core"

    # Check for angular furanocoumarin core
    if mol.HasSubstructMatch(angular_pattern):
        return True, "Contains angular furanocoumarin core"

    # If neither pattern matches, the molecule is not a furanocoumarin
    return False, "Does not contain furanocoumarin core"