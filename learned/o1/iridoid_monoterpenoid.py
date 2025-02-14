"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
"""
Classifies: iridoid monoterpenoid
"""

from rdkit import Chem

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid based on its SMILES string.
    An iridoid monoterpenoid typically consists of a cyclopentane ring fused to a six-membered
    oxygen-containing heterocycle (tetrahydropyran ring), forming a specific bicyclic system.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is an iridoid monoterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the iridoid core SMARTS pattern
    # This pattern represents a cyclopentane ring fused to a tetrahydropyran ring
    iridoid_core_smarts = 'C1CCC2OC1CCC2'
    pattern = Chem.MolFromSmarts(iridoid_core_smarts)

    if pattern is None:
        return False, "Invalid SMARTS pattern for iridoid core"

    # Check if the molecule contains the iridoid core
    if not mol.HasSubstructMatch(pattern):
        return False, "Iridoid core structure not found"

    # Verify monoterpenoid backbone (typically 10 carbon atoms)
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 10:
        return False, f"Too few carbons for monoterpenoid (found {num_carbons}, need at least 10)"

    return True, "Contains iridoid core structure typical of iridoid monoterpenoids"