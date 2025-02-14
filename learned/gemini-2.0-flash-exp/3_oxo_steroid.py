"""
Classifies: CHEBI:47788 3-oxo steroid
"""
from rdkit import Chem

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    A 3-oxo steroid is a steroid with a carbonyl group (=O) at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS for steroid core (more general, focuses on fused rings, including at least 2 six-membered)
    # This pattern identifies fused rings. The [R] means any ring atom, and we want at least one of them to have
    # a 6 membered ring at the start and end, and the intermediate rings can be any size.
    # The pattern must contain at least 3 rings.
    steroid_core_pattern = Chem.MolFromSmarts("[R6]1[R][R]2[R][R]3[R][R]1[R][R]4[R]2[R][R]3[R6]4")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Molecule does not have a steroid core"

    # SMARTS for carbonyl at position 3. It will be attached to a carbon of a 6 membered ring
    # in the steroid core.
    oxo_at_3_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[R6]")
    if not mol.HasSubstructMatch(oxo_at_3_pattern):
        return False, "Molecule does not have a carbonyl at position 3"

    # Verify that the carbonyl is part of the steroid core:
    # Find match of core, then carbonyl, check that the carbon of the carbonyl is in the core
    core_matches = mol.GetSubstructMatches(steroid_core_pattern)
    carbonyl_matches = mol.GetSubstructMatches(oxo_at_3_pattern)
    if not core_matches or not carbonyl_matches:
        return False, "Molecule does not have both a steroid core and a carbonyl"
    
    carbonyl_carbon_in_core = False

    for core_match in core_matches:
        for carbonyl_match in carbonyl_matches:
            carbonyl_carbon_idx = carbonyl_match[0] # the carbon is always the first match
            if carbonyl_carbon_idx in core_match:
                carbonyl_carbon_in_core = True
                break
        if carbonyl_carbon_in_core:
            break
    
    if not carbonyl_carbon_in_core:
        return False, "Carbonyl at 3 not part of the steroid core"
    

    return True, "Molecule is a 3-oxo steroid"