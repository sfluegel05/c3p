"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
from rdkit import Chem

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid based on its SMILES string.
    A short-chain fatty acid is an aliphatic monocarboxylic acid with a chain length of less than C6.
    Any non-hydrocarbon substituent means the compound is not normally regarded as a short-chain fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a short-chain fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for single carboxyl group presence
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if len(carboxyl_matches) != 1:
        return False, "Requires exactly one carboxyl group"

    # Count carbons and ensure no extra elements beyond C, H, and O
    c_count = 0
    extra_elements = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num == 6:  # Carbon
            c_count += 1
        elif atomic_num not in [1, 8]:  # Only Hydrogen and Oxygen are fine
            extra_elements = True

    if extra_elements:
        return False, "Presence of non-hydrocarbon substituents disqualifies compound"

    # Check that the chain length is less than C6
    if c_count >= 6:
        return False, f"Chain length is C{c_count}, should be less than C6"

    # Ensure the structure consists of a single chain (not rings, parallel structures)
    if not mol.GetRingInfo().NumRings() == 0:
        return False, "Rings present, should be an open chain"

    return True, "Valid short-chain fatty acid"

# Example usage
smiles_example = "CCCC(C)C(O)=O"  # 2-methylvaleric acid
result, reason = is_short_chain_fatty_acid(smiles_example)
print(result, reason)