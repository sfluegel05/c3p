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
    Catechins are hydroxyflavans possessing a flavan-3-ol core (i.e. a 2-phenyl-3,4-dihydro-2H-chromene-3-ol skeleton)
    and its substituted derivatives.

    This implementation uses a simplified SMARTS to detect the fused-ring flavan-3-ol core.
    The SMARTS below is defined as: 
        "c1ccc(c(c1))C2C(O)Oc3ccccc23"
    which represents a fused three-ring system in which a phenyl ring (B ring) is attached at position 2, 
    and the chroman ring (A + C rings) carries a hydroxyl at the 3-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a catechin derivative, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string using RDKit.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a simplified SMARTS pattern for a flavan-3-ol core:
    # This pattern looks for a benzopyran fused system with an OH at the appropriate position.
    # Note: Stereochemistry and additional substituents are ignored in this simplified search.
    core_smarts = "c1ccc(c(c1))C2C(O)Oc3ccccc23"
    core_pattern = Chem.MolFromSmarts(core_smarts)
    if core_pattern is None:
        return False, "Failed to create SMARTS pattern for flavan-3-ol core"
    
    # Check if the molecule has a matching substructure for the core.
    if mol.HasSubstructMatch(core_pattern):
        return True, "Molecule contains a flavan-3-ol core (catechin skeleton)"
    else:
        return False, "Molecule does not contain the expected flavan-3-ol core"

# Example testing (uncomment the lines below to run a simple test)
# test_smiles = "O[C@@H]1Cc2c(O)cc(O)cc2O[C@H]1c1ccc(O)c(O)c1"  # (-)-catechin
# result, reason = is_catechin(test_smiles)
# print(result, reason)