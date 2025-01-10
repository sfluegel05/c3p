"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
from rdkit import Chem

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.

    An epoxy fatty acid is characterized by a long carbon chain, an epoxide ring,
    and a carboxylic acid group.

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
    epoxide_pattern = Chem.MolFromSmarts("[C@]1(O[C@H]1)C")
    if not mol.HasSubstructMatch(epoxide_pattern):
        return False, "No epoxide group found"

    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Confirm a chain-like structure appropriate for a fatty acid
    # Could match with 14-30 carbon atoms in linear unsaturated forms
    carbon_chain_pattern = Chem.MolFromSmarts("C(~C(~C(~C(~C~C~C~C~C~C~C~C~C~C~C)=*)=*))*")  # Adjusted for demonstration
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "Insufficient carbon chain length for fatty acid structure"

    # Reassess the number of carbon atoms to ensure a reasonable range for fatty acids
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 14 or c_count > 30:
        return False, f"Unexpected carbon count for a fatty acid, found {c_count}"

    # Use specific functional group and structural detection rather than broad patterns
    complex_undesired = Chem.MolFromSmarts("[R2]")  # Represents complex unwanted ring systems, to be refined
    if mol.HasSubstructMatch(complex_undesired):
        return False, "Complex or polycyclic structures detected"

    return True, "Identified as an epoxy fatty acid with appropriate functional groups and chain length"