"""
Classifies: CHEBI:47787 11-oxo steroid
"""
from rdkit import Chem

def is_11_oxo_steroid(smiles: str):
    """
    Determines if a molecule is an 11-oxo steroid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an 11-oxo steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for the steroid core (simplified).
    # This pattern looks for four fused rings, it's not as specific as the previous
    # one, but will match all steroid cores.
    steroid_core_smarts = "[C]12[C][C]([C]3)[C]4[C]([C]1)[C][C]([C]2)[C]3[C][C]4"
    steroid_core_pattern = Chem.MolFromSmarts(steroid_core_smarts)


    # SMARTS pattern for the 11-oxo group.
    oxo_11_smarts = "[C](-[C])(-[C])-[C](=O)-[C](-[C])(-[C])"
    oxo_11_pattern = Chem.MolFromSmarts(oxo_11_smarts)


    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Not a steroid core"

    # Look for an oxo group at position 11
    if not mol.HasSubstructMatch(oxo_11_pattern):
        return False, "No carbonyl group at position 11"

    return True, "11-oxo steroid"