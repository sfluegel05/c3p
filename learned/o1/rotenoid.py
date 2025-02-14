"""
Classifies: CHEBI:71543 rotenoid
"""
"""
Classifies: rotenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_rotenoid(smiles: str):
    """
    Determines if a molecule is a rotenoid based on its SMILES string.
    A rotenoid contains a cis-fused tetrahydrochromeno[3,4-b]chromene skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a rotenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the rotenoid core structure
    # The pattern represents the cis-fused tetrahydrochromeno[3,4-b]chromene skeleton
    rotenoid_smarts = """
    [#6;R2]1:[#6]:[#6]:[#6]:[#6]:[#6]:1
    -[#6]2-[#6]3-[#6]-[#6]-[#6]-[#6]-3-[#6]-2
    """

    # Remove line breaks and extra whitespace from SMARTS pattern
    rotenoid_smarts = "".join(rotenoid_smarts.strip().split())

    # Create a molecule from the SMARTS pattern
    rotenoid_pattern = Chem.MolFromSmarts(rotenoid_smarts)
    if rotenoid_pattern is None:
        return False, "Invalid SMARTS pattern for rotenoid"

    # Check for substructure match
    if mol.HasSubstructMatch(rotenoid_pattern):
        return True, "Contains the rotenoid core structure"
    else:
        return False, "Does not contain the rotenoid core structure"

__metadata__ = {
    'chemical_class': {
        'name': 'rotenoid',
        'definition': 'Members of the class of tetrahydrochromenochromene that consists of a cis-fused tetrahydrochromeno[3,4-b]chromene skeleton and its substituted derivatives.'
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet'
    },
    'success': True,
    'message': 'Implemented substructure search for rotenoid core'
}