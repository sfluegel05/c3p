"""
Classifies: CHEBI:38757 isoflavones
"""
from rdkit import Chem

def is_isoflavones(smiles: str):
    """
    Determines if a molecule is an isoflavone based on its SMILES string.
    An isoflavone is defined as any isoflavonoid with a 3-aryl-1-benzopyran-4-one
    (3-aryl-4H-chromen-4-one) skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavone, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the isoflavone core SMARTS pattern with variable substitutions
    # The pattern represents the chromen-4-one core with an aryl group at position 3
    isoflavone_smarts = """
    [#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1
    -c2cc(=O)oc3cccc(-c4cccc[#6]4)c23
    """

    # Remove whitespace and newlines from SMARTS pattern
    isoflavone_smarts = isoflavone_smarts.replace('\n', '').replace(' ', '')
    isoflavone_pattern = Chem.MolFromSmarts(isoflavone_smarts)
    if isoflavone_pattern is None:
        return False, "Error in isoflavone core SMARTS pattern"

    # Check if the molecule contains the isoflavone core
    if mol.HasSubstructMatch(isoflavone_pattern):
        return True, "Contains isoflavone core structure"
    else:
        return False, "Does not contain isoflavone core structure"