"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
"""
Classifies: CHEBI:137050 alpha-hydroxy ketone

Definition: A ketone containing a hydroxy group on the alpha-carbon relative to the C=O group.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.
    An alpha-hydroxy ketone is a ketone containing a hydroxy group on the alpha-carbon relative to the carbonyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for alpha-hydroxy ketone
    # [C;H0](=O)[C;H1,H2][O;H1]
    pattern = Chem.MolFromSmarts('[C;H0](=O)[C;H1,H2][O;H1]')

    if mol.HasSubstructMatch(pattern):
        return True, "Contains alpha-hydroxy ketone group"
    else:
        return False, "Does not contain alpha-hydroxy ketone group"