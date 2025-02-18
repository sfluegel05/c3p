"""
Classifies: CHEBI:52221 isothiocyanate
"""
"""
Classifies: CHEBI:52221 isothiocyanate
"""
from rdkit import Chem

def is_isothiocyanate(smiles: str):
    """
    Determines if a molecule is an isothiocyanate based on its SMILES string.
    An isothiocyanate has the general formula R-N=C=S where R is any organic substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isothiocyanate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Enhanced SMARTS pattern with valence checks:
    # [c,C] = aromatic or aliphatic R group
    # [N;D2] = nitrogen with exactly 2 bonds (R and C)
    # [C;D2] = central carbon with exactly 2 bonds (N and S)
    # [S;D1] = sulfur with exactly 1 bond (only to C)
    pattern = Chem.MolFromSmarts(
        "[c,C]-[N;D2]=[C;D2]=[S;D1]"
    )

    # Check for presence of the pattern
    if mol.HasSubstructMatch(pattern):
        return True, "Contains R-N=C=S group with correct valence"
    else:
        return False, "No valid R-N=C=S group detected"