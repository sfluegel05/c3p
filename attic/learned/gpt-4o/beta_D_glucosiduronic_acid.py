"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Glucuronic acid moiety pattern: pyranose ring with hydroxyl groups and carboxylic acid
    glucuronic_acid_pattern = Chem.MolFromSmarts("OC1[C@H]([C@@H](O)[C@H](O)[C@@H](O1)C(=O)O)O")
    if not mol.HasSubstructMatch(glucuronic_acid_pattern):
        return False, "No glucuronic acid moiety found"
    
    # Check for glycosidic bond at the anomeric carbon in beta configuration
    # Searching for beta linkage may involve stereochemistry which can be complex; we will generalize here
    glycosidic_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](OC)O[C@@H](O)[C@H](O)[C@H](O1)C(=O)O")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No beta-glycosidic linkage detected"

    return True, "Contains beta-D-glucosiduronic acid moiety with appropriate glycosidic linkage"