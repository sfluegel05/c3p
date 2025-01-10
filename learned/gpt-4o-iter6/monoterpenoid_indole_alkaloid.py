"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.
    Monoterpenoid indole alkaloids have an indole ring fused with monoterpene-derived structures.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid indole alkaloid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for indole ring system (benzene fused with pyrrole)
    indole_pattern = Chem.MolFromSmarts("c1ccc2[nH]c(cc2c1)")
    if not mol.HasSubstructMatch(indole_pattern):
        return False, "No indole ring found"

    # Check for presence of nitrogen (common in alkaloids)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 1:
        return False, "No nitrogen atom, unlikely to be an alkaloid"

    # Estimate molecular complexity, which is typically high in this class
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 3:
        return False, f"Found {num_rings} rings, not complex enough"

    return True, "Contains indole ring and features typical of monoterpenoid indole alkaloids"