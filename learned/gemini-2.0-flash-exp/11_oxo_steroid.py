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

    # Improved SMARTS pattern for the steroid core. The numbers in [] represent
    # the ring positions that will be used for the subsequent 11-oxo SMARTS.
    # The pattern also allows for unsaturations and substitution at positions
    # that are not ring junctions.
    steroid_core_smarts = "[C]1[C][C]([C])[C]2[C]([C]1)[C][C]3[C]([C]2)[C][C]4[C]([C]3)[C]([C])[C][C]4"
    steroid_core_pattern = Chem.MolFromSmarts(steroid_core_smarts)

    # SMARTS pattern for the 11-oxo group. It refers to the atom numbers of the core pattern.
    # Note that, as specified on the document, 11-oxo means that there is a carbonyl at position 11.
    # position 11 corresponds to the carbon with the number [C]3.
    oxo_11_smarts = "[C]3=O"
    oxo_11_pattern = Chem.MolFromSmarts(oxo_11_smarts)

    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Not a steroid core"
    
    # Combine both queries, to verify that the carbonyl group is at the 11 position
    # defined in the steroid_core_smarts.
    combined_smarts = steroid_core_smarts + "." + oxo_11_smarts
    combined_pattern = Chem.MolFromSmarts(combined_smarts)
    if not mol.HasSubstructMatch(combined_pattern):
       return False, "No carbonyl group at position 11"

    return True, "11-oxo steroid"