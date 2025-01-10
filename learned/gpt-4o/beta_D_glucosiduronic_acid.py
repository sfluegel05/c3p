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

    # Glucuronic acid moiety pattern: a beta-D-glucopyranose with a carboxylic group
    # The pattern should detect the C1 position (anomeric carbon) having beta orientation glycosidic bond
    # Verify both pyranose ring with carboxylic acid expression
    try:
        glucuronic_acid_pattern = Chem.MolFromSmarts("OC1C(O)C(O)C(O)C(O1)C(=O)O")
        if not mol.HasSubstructMatch(glucuronic_acid_pattern):
            return False, "No D-glucuronic acid moiety found"
        
        # Beta configuration of the anomeric carbon in a glycosidic bond
        glycosidic_pattern = Chem.MolFromSmarts("C1(O[C@H]2OC(C(O)=O)C(O)C(O)C2O)C(O)C(O)C(O)C1")
        if not mol.HasSubstructMatch(glycosidic_pattern): 
            return False, "No beta-glycosidic linkage detected"
    except Exception as e:
        return False, f"Error detecting linkage: {str(e)}"

    return True, "Contains beta-D-glucosiduronic acid moiety with appropriate glycosidic linkage"