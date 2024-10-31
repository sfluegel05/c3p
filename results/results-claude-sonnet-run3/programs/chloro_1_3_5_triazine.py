from rdkit import Chem
from rdkit.Chem import AllChem

def is_chloro_1_3_5_triazine(smiles: str):
    """
    Determines if a molecule is a chloro-1,3,5-triazine.
    A chloro-1,3,5-triazine is a 1,3,5-triazine substituted by at least one chloro group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chloro-1,3,5-triazine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Find 1,3,5-triazine core with chloro substituent
    # Using explicit aromatic notation and accounting for different possible orientations
    triazine_chloro_patterns = [
        Chem.MolFromSmarts('n1c(Cl)nc(*)nc1'),
        Chem.MolFromSmarts('n1cnc(Cl)nc1'),
        Chem.MolFromSmarts('n1c(*)nc(Cl)nc1')
    ]
    
    # Check if molecule has triazine core first
    triazine_pattern = Chem.MolFromSmarts('n1cncnc1')
    if not mol.HasSubstructMatch(triazine_pattern):
        return False, "No 1,3,5-triazine core found"
    
    # Check for chloro substituent on triazine
    for pattern in triazine_chloro_patterns:
        if mol.HasSubstructMatch(pattern):
            # Count chloro substituents on triazine
            matches = len(mol.GetSubstructMatches(pattern))
            return True, f"1,3,5-triazine with {matches} chloro substituent(s)"
    
    return False, "1,3,5-triazine found but no chloro substituents"
# Pr=None
# Recall=None