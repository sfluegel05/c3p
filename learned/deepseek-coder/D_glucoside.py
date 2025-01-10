"""
Classifies: CHEBI:35436 D-glucoside
"""
"""
Classifies: CHEBI:16238 D-glucoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.
    A D-glucoside is a molecule where a D-glucose moiety is attached to another group via a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a D-glucoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the D-glucose pattern (beta-D-glucopyranose)
    glucose_pattern = Chem.MolFromSmarts("[C@H]1([C@H]([C@@H]([C@H]([C@@H](O1)CO)O)O)O)")
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No D-glucose moiety found"

    # Check for a glycosidic bond (oxygen connected to the anomeric carbon of glucose)
    glycosidic_bond_pattern = Chem.MolFromSmarts("[C@H]1([C@H]([C@@H]([C@H]([C@@H](O1)CO)O)O)O)-O")
    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No glycosidic bond found"

    # Check that the glucose moiety is not free (i.e., it is attached to another group)
    free_glucose_pattern = Chem.MolFromSmiles("C([C@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)O")
    if mol.HasSubstructMatch(free_glucose_pattern):
        return False, "Glucose moiety is not attached to another group"

    return True, "Contains a D-glucose moiety attached via a glycosidic bond"