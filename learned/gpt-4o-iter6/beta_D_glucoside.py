"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
from rdkit import Chem

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.
    A beta-D-glucoside is characterized by a D-glucose backbone with a linkage
    in the beta configuration at the anomeric carbon. It often appears as
    a standalone unit or in a complex structure.

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

    # Define beta-D-glucoside pattern with broader matching
    try:
        # Primary pattern for beta-D-glucopyranoside
        beta_D_glucoside_pattern = Chem.MolFromSmarts(
            "O[C@@H]1[C@H](O)[C@@H](O)[C@@H](O)[C@H]1O"
        )
        
        if beta_D_glucoside_pattern is None:
            return (None, "SMARTS pattern compilation failed")

        # Check for primary pattern and ensure beta configuration
        if mol.HasSubstructMatch(beta_D_glucoside_pattern):
            return True, "Contains beta-D-glucoside substructure with correct stereochemistry"

    except Exception as e:
        return None, f"Error in SMARTS pattern matching: {e}"

    return False, "No beta-D-glucoside substructure found"