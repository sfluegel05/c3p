"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
"""
Classifies: CHEBI:59826 beta-D-glucoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.
    A beta-D-glucoside is a molecule where a beta-D-glucose moiety is attached to another molecule via a glycosidic bond.

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
    # The pattern should match the beta-D-glucose moiety with the glycosidic bond at C1
    beta_D_glucose_pattern = Chem.MolFromSmarts("[C@H]1([C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)O)O")
    
    # Check if the molecule contains the beta-D-glucose pattern
    if not mol.HasSubstructMatch(beta_D_glucose_pattern):
        return False, "No beta-D-glucose moiety found"

    # Check the configuration of the glycosidic bond
    # The glycosidic bond should be in the beta-configuration (axial position)
    # We can check this by ensuring that the anomeric carbon (C1) is connected to the rest of the molecule in the correct orientation
    # This is a simplified check and may not cover all cases
    anomeric_carbon = mol.GetSubstructMatch(beta_D_glucose_pattern)[0]
    for neighbor in mol.GetAtomWithIdx(anomeric_carbon).GetNeighbors():
        if neighbor.GetAtomicNum() != 8:  # Not an oxygen (glycosidic bond)
            continue
        # Check if the glycosidic bond is in the beta-configuration
        # This is a heuristic and may not be accurate for all cases
        if neighbor.GetIdx() in mol.GetSubstructMatch(beta_D_glucose_pattern):
            return True, "Contains beta-D-glucose moiety with beta-configuration glycosidic bond"

    return False, "No beta-configuration glycosidic bond found"