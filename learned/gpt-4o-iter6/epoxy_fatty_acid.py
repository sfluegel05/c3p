"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is an epoxy fatty acid, False otherwise.
        str: Reason for classification.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the epoxide group pattern
    epoxide_pattern = Chem.MolFromSmarts("[C]1[O][C]1")
    if not mol.HasSubstructMatch(epoxide_pattern):
        return False, "No epoxide group found"

    # Check for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Fatty acids typically have 12-22 carbon atoms
    if c_count < 12 or c_count > 22:
        return False, f"Unexpected carbon count for a fatty acid, found {c_count}"

    # Ensure that the structure is mostly a linear chain with an epoxide
    # Ensure it is not overly cyclic or complex
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings > 1:
        return False, "Unexpected complex ring structure"

    return True, "Identified as an epoxy fatty acid with appropriate functional groups and chain length"