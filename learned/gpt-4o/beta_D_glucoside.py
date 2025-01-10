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

    # SMARTS pattern for a beta-D-glucopyranoside with D-glucose stereochemistry
    # The pattern is refined to specify the beta-anomer at the 1-position (anomeric center) with appropriate stereochemistry
    beta_D_glucoside_pattern = Chem.MolFromSmarts("C1[C@@H](O[C@@H]2[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)O2)[C@@H](O)[C@H](O)[C@H]1O")

    # Check if the molecule has at least one beta-D-glucoside unit
    if not mol.HasSubstructMatch(beta_D_glucoside_pattern):
        return False, "No beta-D-glucose unit with the correct configuration found"

    # There may be additional logic needed to confirm the glycosidic linkage depending on context,
    # For most practical purposes, matching the stereochemistry of the glucose unit constitutes identification.
    
    return True, "Contains beta-D-glucoside with proper beta-configuration at the anomeric center"