"""
Classifies: CHEBI:134251 guaiacols
"""
from rdkit import Chem

def is_guaiacols(smiles: str):
    """
    Determines if a molecule is a guaiacol based on its SMILES string.
    A guaiacol is a phenol with a methoxy group at the ortho position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a guaiacol, False otherwise.
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Create SMARTS pattern for a guaiacol substructure.
    #  'c1' is a ring atom.
    #   [OH] is the phenol
    #   [OCH3] is the methoxy
    #   The ~ are single bonds and the ring bond is implied.
    #   This ensures they are connected at ortho position on the aromatic ring.

    guaiacol_pattern = Chem.MolFromSmarts("c1([OH])c([OCH3])ccccc1") # Improved pattern

    if guaiacol_pattern is None:
       return None, "Invalid SMARTS pattern"

    # Check for the presence of the specific guaiacol substructure.
    if mol.HasSubstructMatch(guaiacol_pattern):
        return True, "Has Guaiacol substructure"
    else:
        return False, "Does not have Guaiacol substructure"