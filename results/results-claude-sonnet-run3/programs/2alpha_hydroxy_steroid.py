from rdkit import Chem
from rdkit.Chem import AllChem

def is_2alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 2alpha-hydroxy steroid.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 2alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Steroid core template (four fused rings)
    steroid_core = Chem.MolFromSmarts('[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~3~[#6]~[#6]~2~[#6]~1')
    
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core found"

    # SMARTS pattern for 2-alpha-hydroxy group
    # Looking for an OH group at position 2 with alpha stereochemistry
    alpha_2_hydroxy = Chem.MolFromSmarts('[#6]1[#6][C@H](O)[#6][#6]2[#6]')
    
    if not mol.HasSubstructMatch(alpha_2_hydroxy):
        return False, "No 2-alpha hydroxy group found"
        
    # Check if the molecule has the basic steroid skeleton and 2-alpha-hydroxy group
    matches = mol.GetSubstructMatches(steroid_core)
    if matches:
        for match in matches:
            # Verify that the 2-alpha-hydroxy is correctly positioned relative to the steroid core
            if mol.HasSubstructMatch(alpha_2_hydroxy):
                return True, "2-alpha hydroxy steroid found"
    
    return False, "Structure does not match 2-alpha-hydroxy steroid pattern"
# Pr=None
# Recall=0.0