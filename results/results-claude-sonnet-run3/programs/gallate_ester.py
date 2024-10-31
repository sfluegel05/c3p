from rdkit import Chem
from rdkit.Chem import AllChem

def is_gallate_ester(smiles: str):
    """
    Determines if a molecule is a gallate ester (ester of gallic acid).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a gallate ester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for galloyl group (3,4,5-trihydroxybenzoyl ester)
    # This pattern matches:
    # - A benzene ring with three hydroxyl groups in 3,4,5 positions
    # - Connected to a carbonyl (C=O)
    # - Which is connected to an oxygen that's connected to something else (ester)
    galloyl_pattern = Chem.MolFromSmarts('[OX2H1]-c1c([OX2H1])c([OX2H1])cc(C(=O)[OX2]-[!H])-c1')
    
    if not mol.HasSubstructMatch(galloyl_pattern):
        return False, "No galloyl ester group found"
    
    matches = mol.GetSubstructMatches(galloyl_pattern)
    num_matches = len(matches)
    
    if num_matches > 0:
        if num_matches == 1:
            return True, "Found 1 galloyl ester group"
        else:
            return True, f"Found {num_matches} galloyl ester groups"
    
    return False, "No valid galloyl ester groups found"
# Pr=None
# Recall=0.0