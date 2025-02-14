"""
Classifies: CHEBI:23824 diol
"""
"""
Classifies: CHEBI:23824 diol
"""
from rdkit import Chem

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol based on its SMILES string.
    A diol is a compound that contains two hydroxy groups, generally assumed to be,
    but not necessarily, alcoholic. Aliphatic diols are also called glycols.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define alcoholic hydroxy group pattern (OH group attached to sp3 carbon, excluding carboxylic acids)
    alcohol_pattern = Chem.MolFromSmarts("[C;!$(C=O)][OX2H]")
    if alcohol_pattern is None:
        return False, "Error in alcohol SMARTS pattern"

    # Find all alcoholic hydroxy groups in the molecule
    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)
    num_alcohol_hydroxy = len(alcohol_matches)

    # Check if the molecule contains at least two alcoholic hydroxy groups
    if num_alcohol_hydroxy >= 2:
        return True, f"Molecule contains {num_alcohol_hydroxy} alcoholic hydroxy groups"
    else:
        return False, f"Molecule contains {num_alcohol_hydroxy} alcoholic hydroxy group(s), less than 2 required for diol"