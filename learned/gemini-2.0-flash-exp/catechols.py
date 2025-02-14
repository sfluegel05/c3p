"""
Classifies: CHEBI:33566 catechols
"""
from rdkit import Chem

def is_catechols(smiles: str):
    """
    Determines if a molecule is a catechol based on its SMILES string.
    A catechol is defined as a compound containing an o-diphenol (two hydroxyl groups on adjacent carbons) component.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catechol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for o-diphenol
    # Pattern: Benzene ring with two adjacent hydroxyl groups
    catechol_pattern = Chem.MolFromSmarts("c1c(O)c(O)cccc1") # or "c1c(O)c(O)cc[cH]1"
    if catechol_pattern is None:
        return None, "Invalid SMARTS pattern"

    # Check for substructure match
    if mol.HasSubstructMatch(catechol_pattern):
        return True, "Molecule contains an ortho-diphenol moiety"
    else:
        return False, "Molecule does not contain an ortho-diphenol moiety"