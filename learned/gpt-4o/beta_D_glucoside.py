"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
from rdkit import Chem

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.
    A beta-D-glucoside is identified by the presence of a beta-D-glucose moiety linked through an ether bond.

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

    # Define SMARTS pattern for a beta-D-glucoside structure
    beta_d_glucoside_pattern = Chem.MolFromSmarts("[C@@H]1(O[C@@H](CO)[C@@H](O)[C@H](O)[C@H]1O)CO")
    
    # Check if the molecule has the beta-D-glucoside pattern
    if not mol.HasSubstructMatch(beta_d_glucoside_pattern):
        return False, "No beta-D-glucoside substructure found"

    return True, "Contains beta-D-glucoside substructure"

# Example usage:
# smiles_example = "OC[C@H]1O[C@@H](OCC=2C=CC=CC2)[C@@H]([C@H]([C@@H]1O)O)O"
# result, reason = is_beta_D_glucoside(smiles_example)
# print(result, reason)