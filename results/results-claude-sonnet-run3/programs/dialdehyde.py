from rdkit import Chem
from rdkit.Chem import AllChem

def is_dialdehyde(smiles: str):
    """
    Determines if a molecule is a dialdehyde (contains exactly 2 aldehyde groups).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a dialdehyde, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
        
    # Count aldehyde groups
    # Aldehyde pattern: carbon with double bonded oxygen and single bonded hydrogen
    patt = Chem.MolFromSmarts('[CH](=O)')
    matches = mol.GetSubstructMatches(patt)
    
    num_aldehydes = len(matches)
    
    if num_aldehydes == 0:
        return False, "No aldehyde groups found"
    elif num_aldehydes == 1:
        return False, "Only one aldehyde group found"
    elif num_aldehydes == 2:
        return True, "Two aldehyde groups found"
    else:
        return False, f"More than two aldehyde groups found ({num_aldehydes})"
# Pr=1.0
# Recall=1.0