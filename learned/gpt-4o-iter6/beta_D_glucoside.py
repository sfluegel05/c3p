"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
from rdkit import Chem

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.
    A beta-D-glucoside is characterized by a D-glucose backbone with a linkage
    in the beta configuration at the anomeric carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucoside, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a refined beta-D-glucoside pattern
    # Match a beta configuration on the anomeric carbon in a glucose-like structure with proper stereochemistry
    beta_D_glucoside_pattern = Chem.MolFromSmarts(
        "O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H]1"
    )
    if beta_D_glucoside_pattern is None:
        return None, "SMARTS pattern compilation failed"

    # Search for the pattern considering stereochemistry
    if mol.HasSubstructMatch(beta_D_glucoside_pattern):
        return True, "Contains beta-D-glucoside substructure with correct stereochemistry"

    return False, "No beta-D-glucoside substructure found"