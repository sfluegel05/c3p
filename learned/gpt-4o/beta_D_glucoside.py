"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
from rdkit import Chem

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.
    A beta-D-glucoside is characterized by a beta-D-glucose moiety linked through an ether bond.
    The stereochemistry must reflect the beta-configuration at the anomeric site.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a beta-D-glucoside, False otherwise.
        str: Reason for classification.
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for beta-D-glucoside structures
    # Capture the 6-membered pyranose ring with correct stereochemistry
    # and an ether linkage to an external group at the anomeric carbon.
    beta_d_glucoside_patterns = [
        Chem.MolFromSmarts("[C@H]1(O[C@H](CO)O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)"),
        Chem.MolFromSmarts("O[C@H]1[C@H](O[C@@H](O)[C@@H](O)[C@H](O)[C@H]1O)CO"),
        Chem.MolFromSmarts("O[C@H]1[C@@H](O[C@H](O)[C@@H](O)[C@H](O)[C@H]1O)CO")
    ]

    # Check if the molecule contains the beta-D-glucoside substructure
    for pattern in beta_d_glucoside_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains beta-D-glucoside substructure"
    
    return False, "No beta-D-glucoside substructure found"

# Example usage:
# result, reason = is_beta_D_glucoside("OC[C@H]1O[C@@H](Oc2ccc3ccc(c2)-c2ccccc3c2)C1O")
# print(result, reason)