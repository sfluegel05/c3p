"""
Classifies: CHEBI:32863 secondary amine
"""
"""
Classifies: secondary amine
"""
from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule is a secondary amine based on its SMILES string.
    A secondary amine is a nitrogen atom bonded to exactly two carbon atoms and one hydrogen,
    excluding cases where nitrogen is part of aromatic systems, amides, or other functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary amine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for secondary amine
    # Nitrogen bonded to one hydrogen and two carbons (non-carbonyl), not aromatic
    secondary_amine_smarts = "[#7H1;!a;$([N;H1][C;!$(C=O)][C;!$(C=O)])]"

    pattern = Chem.MolFromSmarts(secondary_amine_smarts)
    if pattern is None:
        return False, "Invalid SMARTS pattern"

    matches = mol.GetSubstructMatches(pattern)
    if matches:
        return True, "Contains secondary amine group"
    else:
        return False, "Does not contain secondary amine group"