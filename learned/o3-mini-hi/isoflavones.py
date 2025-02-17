"""
Classifies: CHEBI:38757 isoflavones
"""
"""
Classifies: isoflavones – Any isoflavonoid with a 3-aryl-1-benzopyran-4-one (3-aryl-4H-chromen-4-one) skeleton and its substituted derivatives.
"""

from rdkit import Chem

def is_isoflavones(smiles: str):
    """
    Determines if a molecule belongs to the isoflavones group based on its SMILES string.
    Isoflavones are defined as isoflavonoids containing a 3-aryl-1-benzopyran-4-one core.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule has an isoflavone core, False otherwise
        str: Reason for classification or rejection
    """
    # Parse the SMILES string using RDKit.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # We define a SMARTS pattern that aims to capture the 3-aryl-1-benzopyran-4-one (isoflavone) core.
    # Explanation:
    #   The pattern "c1ccc2c(c1)oc(=O)c(c2)-c3ccccc3" breaks down as follows:
    #     • "c1ccc2c(c1)oc(=O)c(c2)" represents the fused benzopyran-4-one skeleton.
    #     • "-c3ccccc3" requires that a phenyl ring (the aryl substituent) is attached 
    #       to one of the carbons in the pyran ring (here, representing the C-3 position).
    # This pattern should be present in molecules such as formononetin, glycitein, genistein, etc.
    isoflavone_smarts = "c1ccc2c(c1)oc(=O)c(c2)-c3ccccc3"
    core = Chem.MolFromSmarts(isoflavone_smarts)
    if core is None:
        return False, "Failed to generate SMARTS pattern for isoflavone core"

    # Check if the molecule contains this isoflavone substructure.
    if not mol.HasSubstructMatch(core):
        return False, "Molecule does not contain the 3-aryl-1-benzopyran-4-one (isoflavone) core"

    return True, "Molecule contains the 3-aryl-1-benzopyran-4-one (isoflavone) core structure"