"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

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
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "Requires exactly one carboxyl group"
    
    # Ensure chain is open (no rings)
    if mol.GetRingInfo().NumRings() != 0:
        return False, "Rings present, should be an open chain"
    
    # Count non-hydrocarbon elements (other than C, H, O)
    num_others = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in {1, 6, 8})
    if num_others > 0:
        return False, "Presence of non-hydrocarbon substituents disqualifies compound"

    # Compute the longest carbon chain
    carbon_atom_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    carbon_paths = [rdmolops.GetShortestPath(mol, carbon_atom_indices[0], carbon_atom_indices[-1])]
    longest_chain_length = max(len(path) for path in carbon_paths)

    if longest_chain_length <= 6:
        return True, "Valid short-chain fatty acid"
    else:
        return False, f"Carbon chain length is {longest_chain_length}, should be less than C6"

# Example usage
smiles_example = "CCCC(O)=O"  # Butyric acid
result, reason = is_short_chain_fatty_acid(smiles_example)
print(result, reason)