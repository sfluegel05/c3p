"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
from rdkit import Chem

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.
    A beta-D-glucoside is characterized by a beta-configuration at the anomeric carbon of a D-glucose unit.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains a beta-D-glucoside, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into an RDKit Mol object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern to match the beta-D-glucopyranoside structure. 
    beta_D_glucopyranoside_pattern = Chem.MolFromSmarts("O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)[C@@H]1O")

    # Check if the molecule has at least one beta-D-glucoside unit.
    if not mol.HasSubstructMatch(beta_D_glucopyranoside_pattern):
        return False, "No beta-D-glucose unit with the correct configuration found"

    return True, "Contains beta-D-glucoside with proper beta-configuration at the anomeric center"