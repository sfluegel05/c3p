from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_apo_carotenoid_triterpenoid(smiles: str):
    """
    Determines if a molecule is an apo carotenoid triterpenoid.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an apo carotenoid triterpenoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
        
    # Check for minimum number of carbons (typically >20)
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    if num_carbons < 20:
        return False, "Too few carbons for an apo carotenoid triterpenoid"
        
    # Check for conjugated double bond system
    conjugated_db = False
    matches = mol.GetSubstructMatches(Chem.MolFromSmarts('C=CC=C'))
    if len(matches) > 2:
        conjugated_db = True
        
    if not conjugated_db:
        return False, "Lacks characteristic conjugated double bond system"
        
    # Check for characteristic terminal groups
    terminal_groups = False
    
    # Check for carboxyl group
    if mol.HasSubstructMatch(Chem.MolFromSmarts('C(=O)O')):
        terminal_groups = True
        
    # Check for aldehyde group
    if mol.HasSubstructMatch(Chem.MolFromSmarts('C=O')):
        terminal_groups = True
        
    # Check for alcohol group
    if mol.HasSubstructMatch(Chem.MolFromSmarts('CO')):
        terminal_groups = True
        
    if not terminal_groups:
        return False, "Missing characteristic terminal groups"
        
    # Check for methyl branching pattern
    methyl_branches = len(mol.GetSubstructMatches(Chem.MolFromSmarts('CC(C)=C')))
    if methyl_branches < 2:
        return False, "Insufficient methyl branching pattern"
        
    return True, "Contains characteristic conjugated system, terminal groups and methyl branching pattern typical of apo carotenoids"
# Pr=1.0
# Recall=0.7142857142857143