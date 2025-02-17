"""
Classifies: CHEBI:23053 catechin
"""
"""
Classifies: Catechins (Members of the class of hydroxyflavan that have a flavan-3-ol skeleton
and its substituted derivatives).
"""

from rdkit import Chem

def is_catechin(smiles: str):
    """
    Determines if a molecule belongs to the catechin class based on its SMILES string.
    Catechins are hydroxyflavans possessing a flavan-3-ol core (a 2-phenyl-3,4-dihydro-2H-chromene-3-ol skeleton)
    and its substituted derivatives.
    
    To improve the detection (and allow for substituted derivatives or slight variations in saturation),
    we first remove all stereochemistry and then use a relaxed SMARTS for the catechin core.
    The relaxed SMARTS is defined as:
        "c1ccc(c(c1))C2CC(O)c3ccccc3O2"
    which represents a B ring (phenyl) attached to a chroman-like system carrying an –OH at the appropriate position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if visible catechin core is detected, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Remove stereochemistry to simplify the matching process.
    Chem.RemoveStereochemistry(mol)
    
    # Define a relaxed SMARTS pattern for the flavan-3-ol (catechin) core.
    # This pattern handles a 2-phenyl-3,4-dihydro-2H-chromene-3-ol structure.
    core_smarts = "c1ccc(c(c1))C2CC(O)c3ccccc3O2"
    core_pattern = Chem.MolFromSmarts(core_smarts)
    if core_pattern is None:
        return False, "Failed to create SMARTS pattern for catechin core"

    # Remove stereochemistry from the core pattern also.
    Chem.RemoveStereochemistry(core_pattern)
    
    # Check for substructure match to detect the catechin core.
    if mol.HasSubstructMatch(core_pattern):
        return True, "Molecule contains a catechin flavan-3-ol core (with potential substituents)"
    else:
        return False, "Molecule does not contain the expected flavan-3-ol core"

# Example testing (uncomment below lines to run a simple test)
# test_smiles = "O[C@@H]1Cc2c(O)cc(O)cc2O[C@H]1c1ccc(O)c(O)c1"  # (-)-catechin
# result, reason = is_catechin(test_smiles)
# print(result, reason)