"""
Classifies: CHEBI:23824 diol
"""
"""
Classifies: CHEBI:23824 diol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Define hydroxy group pattern (Oxygen with single bonds and one attached hydrogen)
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    if hydroxy_pattern is None:
        return False, "Error in hydroxy SMARTS pattern"

    # Find all hydroxy groups in the molecule
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    num_hydroxy = len(hydroxy_matches)

    # Check if the molecule contains at least two hydroxy groups
    if num_hydroxy >= 2:
        return True, f"Molecule contains {num_hydroxy} hydroxy groups"
    else:
        return False, f"Molecule contains {num_hydroxy} hydroxy group(s), less than required for diol"