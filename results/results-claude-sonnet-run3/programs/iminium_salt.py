from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import IPythonConsole

def is_iminium_salt(smiles: str):
    """
    Determines if a molecule is an iminium salt (R2C=N(+)R2).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an iminium salt, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
        
    # Look for positively charged nitrogen atoms
    iminium_patterns = []
    
    # Pattern 1: Basic iminium pattern R2C=N+(R)R
    pattern1 = Chem.MolFromSmarts('[C]=[N+;!$(N-[O,N,S])]')
    if pattern1 is not None and mol.HasSubstructMatch(pattern1):
        iminium_patterns.append("R2C=N+(R)R")
        
    # Pattern 2: Iminium pattern in aromatic ring
    pattern2 = Chem.MolFromSmarts('[c]:[n+;!$(n-[O,N,S])]')
    if pattern2 is not None and mol.HasSubstructMatch(pattern2):
        iminium_patterns.append("aromatic C=N+")
        
    # Check for counterions
    anion_patterns = [
        Chem.MolFromSmarts('[Cl-,Br-,I-,BF4-,PF6-]'),
        Chem.MolFromSmarts('[Fe-]')
    ]
    
    has_counterion = False
    for pattern in anion_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            has_counterion = True
            break
            
    if len(iminium_patterns) > 0 and has_counterion:
        return True, f"Iminium salt with pattern(s): {', '.join(iminium_patterns)}"
    elif len(iminium_patterns) > 0 and not has_counterion:
        return False, "Has iminium cation but no counterion found"
    else:
        return False, "No iminium cation pattern found"
# Pr=1.0
# Recall=0.6666666666666666