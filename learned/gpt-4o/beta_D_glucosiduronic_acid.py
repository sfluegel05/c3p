"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
from rdkit import Chem

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid derivative based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucosiduronic acid derivative, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Glucuronic acid moiety pattern: a D-pyranose ring with hydroxyl groups and carboxylic acid
    glucuronic_acid_pattern = Chem.MolFromSmarts("OC1[C@H]([C@@H](O)[C@H](O)[C@@H](O1)C(=O)O)")
    if not mol.HasSubstructMatch(glucuronic_acid_pattern):
        return False, "No D-glucuronic acid moiety found"
    
    # Beta configuration glycosidic bond involving the anomeric carbon
    try:
        glycosidic_pattern = Chem.MolFromSmarts("[C@H]1(O[C@H])(O)[C@H](O)[C@H](O)[C@@H](O1)C(=O)O")
        if not mol.HasSubstructMatch(glycosidic_pattern):
            return False, "No beta-glycosidic linkage detected"
    except Exception as e:
        return False, f"Error detecting linkage: {str(e)}"

    return True, "Contains beta-D-glucosiduronic acid moiety with appropriate glycosidic linkage"