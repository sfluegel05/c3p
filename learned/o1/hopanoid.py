"""
Classifies: CHEBI:51963 hopanoid
"""
"""
Classifies: CHEBI:27439 hopanoid
"""
from rdkit import Chem

def is_hopanoid(smiles: str):
    """
    Determines if a molecule is a hopanoid based on its SMILES string.
    A hopanoid is a triterpenoid based on a hopane skeleton (pentacyclic triterpene).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hopanoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the hopane skeleton SMARTS pattern (pentacyclic ring system)
    # The pattern represents five fused rings: four cyclohexane and one cyclopentane
    hopane_smarts = "[#6]1[#6][#6][#6]2[#6]([#6]1)[#6][#6][#6]3[#6]([#6]2)[#6][#6][#6]4[#6]([#6]3)[#6][#6][#6]5[#6]([#6]4)[#6][#6][#6]5"
    hopane_mol = Chem.MolFromSmarts(hopane_smarts)
    if hopane_mol is None:
        return False, "Invalid hopane skeleton SMARTS"

    # Use substructure matching without considering stereochemistry
    if not mol.HasSubstructMatch(hopane_mol):
        return False, "No hopane skeleton found"

    return True, "Contains hopane skeleton characteristic of hopanoids"