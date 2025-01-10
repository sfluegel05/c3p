"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
from rdkit import Chem

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid based on its SMILES string.
    An iridoid monoterpenoid typically contains a cyclopentanopyran framework.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule matches the iridoid monoterpenoid class, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved SMARTS pattern for cyclopentanopyran systems, with variations to capture more examples
    # This is generalized to capture variations in the cyclopentanopyran core structure
    cyclopentanopyran_smarts_variants = [
        "C1CCC2C(C1)OC=C2",  # Basic cyclic framework
        "[C&R1]1[C&R1][C&R1][C&R1][C&R1]2C(C1)[O&R1][C&R1]=C2",  # More generalized variant
        "[C@@R1&X4]1[C@R1&X4][C@@R1&X4][C@R1&X4][C@@R1&X4]2[C@R1&X4]([C@R1&X4]1)[O&R1][C&R1]=C2"  # Including stereochemistry
    ]

    # Attempt to match any of the SMARTS patterns to the molecule
    for pattern in cyclopentanopyran_smarts_variants:
        cyclopentanopyran_pattern = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(cyclopentanopyran_pattern):
            return True, "Compound contains a cyclopentanopyran framework indicative of iridoid monoterpenoids"
    
    return False, "No cyclopentanopyran framework found; does not match typical iridoid monoterpenoid structure"