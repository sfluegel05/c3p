from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_thromboxane(smiles: str):
    """
    Determines if a molecule is a thromboxane based on structural features.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a thromboxane, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
        
    # Check for presence of oxane ring (6-membered oxygen heterocycle)
    oxane_pattern = Chem.MolFromSmarts('C1OCCCC1')
    if not mol.HasSubstructMatch(oxane_pattern):
        return False, "No oxane ring found"
        
    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts('C(=O)[OH]')
    
    # Check for hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts('[OH]')
    
    # Check for alkene groups
    alkene_pattern = Chem.MolFromSmarts('C=C')
    
    # Basic thromboxane features include:
    # - Oxane ring
    # - At least one hydroxyl group
    # - Alkene functionality
    # - Often (but not always) a carboxylic acid group
    
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl groups found"
        
    if not mol.HasSubstructMatch(alkene_pattern):
        return False, "No alkene groups found"
        
    # Additional structural features that help identify thromboxanes
    bicyclic_pattern = Chem.MolFromSmarts('C1OC2CC1O2') # For TXA-type
    monocyclic_pattern = Chem.MolFromSmarts('OC1CC(O)C([*])C([*])O1') # For TXB-type
    
    if mol.HasSubstructMatch(bicyclic_pattern):
        return True, "Thromboxane A-type structure identified"
    elif mol.HasSubstructMatch(monocyclic_pattern):
        return True, "Thromboxane B-type structure identified"
    elif mol.HasSubstructMatch(carboxylic_pattern):
        # If it has the basic features plus carboxylic acid, likely a thromboxane derivative
        return True, "Thromboxane derivative with carboxylic acid group"
    else:
        # If it has oxane ring, hydroxyls and alkenes but doesn't match known patterns,
        # it could be a modified thromboxane
        return True, "Possible modified thromboxane structure"
# Pr=0.9375
# Recall=0.8823529411764706