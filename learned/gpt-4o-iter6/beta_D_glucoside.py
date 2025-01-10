"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.
    A beta-D-glucoside is a D-glucose molecule with a beta glycosidic linkage.

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

    # Define a beta-D-glucoside pattern
    # This pattern matches a glucose ring with beta anomeric oxygen linkage
    beta_D_glucoside_pattern = Chem.MolFromSmarts("O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1")
    if beta_D_glucoside_pattern is None:
        return None, "SMARTS pattern compilation failed"

    # Look for beta-D-glucoside pattern
    if not mol.HasSubstructMatch(beta_D_glucoside_pattern):
        return False, "No beta-D-glucoside substructure found"

    return True, "Contains beta-D-glucoside substructure"