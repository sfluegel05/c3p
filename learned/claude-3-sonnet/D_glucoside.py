"""
Classifies: CHEBI:35436 D-glucoside
"""
"""
Classifies: D-glucoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.
    A D-glucoside must have a D-glucose moiety connected to another molecule via a glycosidic bond.

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

    # SMARTS patterns for D-glucose in pyranose form (both alpha and beta)
    # The pattern matches a 6-membered ring with 4 OH groups and 1 CH2OH group
    # The stereochemistry is specified to match D-glucose specifically
    alpha_glucose = Chem.MolFromSmarts('[C@@H]1([OH])O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O')
    beta_glucose = Chem.MolFromSmarts('[C@@H]1([OH])O[C@@H](CO)[C@@H](O)[C@H](O)[C@H]1O')
    
    # Look for glycosidic linkage patterns
    o_glycoside = Chem.MolFromSmarts('[C@@H]1([OR2])O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O')  # O-glycoside
    n_glycoside = Chem.MolFromSmarts('[C@@H]1([NR2])O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O')  # N-glycoside
    s_glycoside = Chem.MolFromSmarts('[C@@H]1([SR2])O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O')  # S-glycoside
    c_glycoside = Chem.MolFromSmarts('[C@@H]1([CR2])O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O')  # C-glycoside

    # Check if molecule contains D-glucose
    has_alpha = mol.HasSubstructMatch(alpha_glucose) if alpha_glucose else False
    has_beta = mol.HasSubstructMatch(beta_glucose) if beta_glucose else False
    
    if not (has_alpha or has_beta):
        return False, "No D-glucose moiety found"

    # Check for glycosidic linkages
    has_o_glycoside = mol.HasSubstructMatch(o_glycoside) if o_glycoside else False
    has_n_glycoside = mol.HasSubstructMatch(n_glycoside) if n_glycoside else False
    has_s_glycoside = mol.HasSubstructMatch(s_glycoside) if s_glycoside else False
    has_c_glycoside = mol.HasSubstructMatch(c_glycoside) if c_glycoside else False

    if not (has_o_glycoside or has_n_glycoside or has_s_glycoside or has_c_glycoside):
        return False, "No glycosidic linkage found"

    # Determine the type of glycoside
    glycoside_types = []
    if has_o_glycoside:
        glycoside_types.append("O-glycoside")
    if has_n_glycoside:
        glycoside_types.append("N-glycoside")
    if has_s_glycoside:
        glycoside_types.append("S-glycoside")
    if has_c_glycoside:
        glycoside_types.append("C-glycoside")

    anomeric_form = "alpha" if has_alpha else "beta"
    
    return True, f"Contains {anomeric_form}-D-glucose in {', '.join(glycoside_types)} form"