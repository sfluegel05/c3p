from rdkit import Chem
from rdkit.Chem import AllChem

def is_methoxy_1_3_5_triazine(smiles: str):
    """
    Determines if a molecule is a methoxy-1,3,5-triazine (1,3,5-triazine with at least one methoxy substituent).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methoxy-1,3,5-triazine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # First check for 1,3,5-triazine core
    triazine_pattern = Chem.MolFromSmarts('[n]1[c][n][c][n][c]1')  # 1,3,5-triazine pattern
    if not mol.HasSubstructMatch(triazine_pattern):
        return False, "No 1,3,5-triazine core found"

    # Check for methoxy group attached to triazine
    methoxy_on_triazine = Chem.MolFromSmarts('[n]1[c]([O][C])[n][c][n][c]1')
    if not mol.HasSubstructMatch(methoxy_on_triazine):
        return False, "No methoxy group attached to triazine ring"

    # Count number of methoxy groups on triazine
    matches = mol.GetSubstructMatches(methoxy_on_triazine)
    num_methoxy = len(matches)

    return True, f"1,3,5-triazine with {num_methoxy} methoxy group(s)"
# Pr=None
# Recall=None