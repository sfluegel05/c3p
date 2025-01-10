"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
from rdkit import Chem

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

    # Look for the epoxide group pattern (three-membered cyclic ether: C1OC1)
    epoxide_pattern = Chem.MolFromSmarts("[C]1-[O]-[C]1")
    if not mol.HasSubstructMatch(epoxide_pattern):
        return False, "No epoxide group found"

    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Confirm a chain-like structure appropriate for a fatty acid
    carbon_chain_pattern = Chem.MolFromSmarts("C~C~C~C~C~C")  # General long carbon chain pattern
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "Insufficient carbon chain length for fatty acid structure"

    # Reassess the number of carbon atoms to ensure a reasonable range for fatty acids
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 14 or c_count > 30:  # Adjust heuristic for upper and lower limits
        return False, f"Unexpected carbon count for a fatty acid, found {c_count}"

    # Avoid configurations indicating ring systems beyond simple alkyl chains
    unwanted_ring_pattern = Chem.MolFromSmarts("C1-=C-=C1")  # Example, complex or polycyclic structures
    if mol.HasSubstructMatch(unwanted_ring_pattern):
        return False, "Unwanted ring or fused ring structures detected"

    return True, "Identified as an epoxy fatty acid with appropriate functional groups and chain length"