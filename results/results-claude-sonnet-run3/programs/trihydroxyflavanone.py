from rdkit import Chem
from rdkit.Chem import AllChem

def is_trihydroxyflavanone(smiles: str):
    """
    Determines if a molecule is a trihydroxyflavanone (a flavanone with exactly 3 hydroxyl groups)
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a trihydroxyflavanone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
        
    # Check for flavanone core structure
    # SMARTS pattern for flavanone core (without hydroxyls)
    flavanone_pattern = Chem.MolFromSmarts('O1c2ccccc2C(=O)CC1c1ccccc1')
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "Not a flavanone core structure"
        
    # Count hydroxyl groups
    oh_pattern = Chem.MolFromSmarts('[OH]')
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    
    num_oh = len(oh_matches)
    
    if num_oh != 3:
        return False, f"Has {num_oh} hydroxyl groups, needs exactly 3"
        
    return True, "Has flavanone core and exactly 3 hydroxyl groups"
# Pr=1.0
# Recall=0.6666666666666666