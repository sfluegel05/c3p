"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
"""
Classifies: CHEBI:59805 beta-D-glucoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.
    A beta-D-glucoside is a molecule where a beta-D-glucose is attached to another molecule via a glycosidic bond.

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

    # Define the beta-D-glucose pattern
    beta_D_glucose_pattern = Chem.MolFromSmarts("[C@@H]1([C@H]([C@@H]([C@H]([C@@H](O1)CO)O)O)O)O")
    if not mol.HasSubstructMatch(beta_D_glucose_pattern):
        return False, "No beta-D-glucose moiety found"

    # Check if the glucose is attached to another molecule via a glycosidic bond
    glycosidic_bond_pattern = Chem.MolFromSmarts("[C@@H]1([C@H]([C@@H]([C@H]([C@@H](O1)CO)O)O)O)O-[*]")
    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No glycosidic bond found"

    # Verify the beta configuration by checking the anomeric carbon
    anomeric_carbon_pattern = Chem.MolFromSmarts("[C@@H]1([C@H]([C@@H]([C@H]([C@@H](O1)CO)O)O)O)O")
    matches = mol.GetSubstructMatches(anomeric_carbon_pattern)
    for match in matches:
        anomeric_carbon = mol.GetAtomWithIdx(match[0])
        if anomeric_carbon.GetHybridization() != Chem.HybridizationType.SP3:
            return False, "Anomeric carbon is not in the correct configuration"

    return True, "Contains a beta-D-glucose moiety attached via a glycosidic bond"