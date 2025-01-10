"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, AllChem

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.
    A beta-D-glucoside features a D-glucose unit with the anomeric carbon in the beta configuration.

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
    
    # SMARTS pattern for beta-D-glucopyranose
    beta_D_glucose_pattern = Chem.MolFromSmarts("O[C@@H]1[C@@H]([C@H]([C@@H]([C@H]([C@H]1O)CO)O)O)O")
    
    if not mol.HasSubstructMatch(beta_D_glucose_pattern):
        return False, "No beta-D-glucose unit with the correct configuration found"
    
    # Check for glycosidic bond (O-linked to another part of the molecule)
    # Assuming the sugar oxygen connects to an anomeric carbon
    if not mol.HasSubstructMatch(Chem.MolFromSmarts("O[C@H]1[C@H](O)C[C@@H](O)C[C@@H]1O")):
        return False, "No glycosidic bond found connecting the beta-D-glucose"
    
    return True, "Contains beta-D-glucose with a glycosidic linkage"