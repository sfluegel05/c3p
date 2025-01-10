"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # SMARTS pattern to detect a beta-D-glucopyranoside
    beta_D_glucoside_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1")

    # Check if the molecule has at least one beta-D-glucopyranoside unit
    if not mol.HasSubstructMatch(beta_D_glucoside_pattern):
        return False, "No beta-D-glucose unit with the correct configuration found"

    # Check for the presence of a glycosidic bond involving the anomeric oxygen (assuming it is linked to another moiety)
    glycosidic_bond_pattern = Chem.MolFromSmarts("O[C@@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1[*]")
    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No glycosidic bond found involving the beta-D-glucose unit"

    return True, "Contains beta-D-glucoside with a glycosidic linkage"