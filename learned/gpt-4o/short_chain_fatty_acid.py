"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
from rdkit import Chem

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid based on its SMILES string.
    A short-chain fatty acid is defined as an aliphatic monocarboxylic acid with a 
    total carbon count in the main chain ≤ 6, with no additional functional groups
    outside the carboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a short-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for exactly one carboxylic acid group: -C(=O)O
    carboxylic_pattern = Chem.MolFromSmarts('C(=O)[O;H1]')
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_pattern)
    if len(carboxylic_matches) != 1:
        return False, f"Contains {len(carboxylic_matches)} carboxylic acid groups, expected exactly 1"

    # Count carbon atoms in the main chain (considering possible branching)
    c_count = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if c_count > 6:
        return False, f"Contains {c_count} carbons, expected 6 or fewer"
    
    # Ensure no non-hydrocarbon substituents (excluding carboxylic group)
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in {'H', 'C', 'O'}:
            return False, f"Contains non-hydrocarbon substituent: {atom.GetSymbol()}"
    
    # Ensure there are no aromatic structures, i.e., aliphatic only
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic structures, expected aliphatic only"
    
    return True, "Valid short-chain fatty acid with carboxylic acid group and ≤6 carbon atoms"

# Example call 
# This should return (True, "Valid short-chain fatty acid with carboxylic acid group and ≤6 carbon atoms")
# print(is_short_chain_fatty_acid("CCCC(O)=O"))  # Represents butyric acid