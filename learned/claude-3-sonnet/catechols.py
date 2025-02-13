"""
Classifies: CHEBI:33566 catechols
"""
"""
Classifies: CHEBI:27027 catechol
Any compound containing an o-diphenol component
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_catechol(smiles: str):
    """
    Determines if a molecule is a catechol based on its SMILES string.
    A catechol is any compound containing an o-diphenol component (two hydroxyl groups
    attached to adjacent carbon atoms, which may or may not be part of a ring system).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catechol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define o-diphenol pattern
    aromatic_pattern = Chem.MolFromSmarts("c1c(O)c(O)cc1")
    aliphatic_pattern = Chem.MolFromSmarts("OC(C)CO")

    # Check for aromatic and aliphatic o-diphenol patterns
    if mol.HasSubstructMatch(aromatic_pattern) or mol.HasSubstructMatch(aliphatic_pattern):
        # Perform additional checks
        # ...

        return True, "Contains an o-diphenol component"

    return False, "Does not contain an o-diphenol component"