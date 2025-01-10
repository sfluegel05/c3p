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

    # Define beta-D-glucoside pattern, considering anomeric C and stereochemistry
    try:
        # Primary pattern for beta-D-glucopyranoside fragment (furanoside not considered here)
        # The beta linkage semantics generally imply "Up" relative stereochemistry for C1-O
        beta_D_glucoside_pattern = Chem.MolFromSmarts(
            "[C@H]1([C@@H](O)[C@H](O)[C@@H](CO)[C@H](O1)O)O"
        )
        
        if beta_D_glucoside_pattern is None:
            return (None, "SMARTS pattern compilation failed")

        # Check for the pattern and ensure beta configuration
        if mol.HasSubstructMatch(beta_D_glucoside_pattern):
            return True, "Contains beta-D-glucoside substructure with correct stereochemistry"

    except Exception as e:
        return None, f"Error in SMARTS pattern matching: {e}"

    return False, "No beta-D-glucoside substructure found"