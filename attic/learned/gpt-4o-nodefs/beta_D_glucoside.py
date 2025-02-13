"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
from rdkit import Chem

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.
    
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

    # Define the SMARTS pattern for the beta-D-glucoside moiety
    # This captures beta linkage at the anomeric position
    beta_D_glucose_pattern = Chem.MolFromSmarts("O[C@@H]1[C@H](O)[C@H](O)[C@H](O)[C@H](CO)O1")

    # Check for the presence of the beta-D-glucoside pattern
    if not mol.HasSubstructMatch(beta_D_glucose_pattern):
        return False, "No beta-D-glucoside moiety found"
    
    # If pattern is found, return True
    return True, "Beta-D-glucoside moiety is present"

# Example usage
smiles_example = "COc1cc2OC[C@H]3Oc4c5C[C@@H](Oc5ccc4[C@H](O)[C@H]3c2cc1OC)C(=C)CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O"
is_beta_D_glucoside(smiles_example)