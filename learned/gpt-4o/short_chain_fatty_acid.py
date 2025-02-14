"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
from rdkit import Chem

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid based on its SMILES string.
    A short-chain fatty acid is defined as an aliphatic monocarboxylic acid with a
    total carbon count excluding the carboxylic group ≤ 6, allowing some functional groups.

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
    
    # Check the main carbon chain length
    # Consider atoms excluding the carboxylic carbon (which is typically seen separately in aliphatic R-COOH)
    non_oxo_carbons = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            if not atom.IsInRing():
                non_oxo_carbons += 1
                for nbr in atom.GetNeighbors():
                    if nbr.GetAtomicNum() != 6 and nbr.GetAtomicNum() != 8:
                        return False, f"Contains non-allowed substituent: {nbr.GetSymbol()}"
    
    # Validate that the carbon backbone contains ≤ 5 carbons if subtracting the carboxy carbon
    if non_oxo_carbons - 1 > 5:
        return False, f"Main chain exceeds C5, found {non_oxo_carbons - 1} carbons excluding carboxylic carbon"
    
    # Ensure no aromatic structures
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic structures, expected aliphatic only"
    
    return True, "Valid short-chain fatty acid with carboxylic acid group and ≤6 carbon atoms"

# Example call 
# This should return (True, "Valid short-chain fatty acid with carboxylic acid group and ≤6 carbon atoms")
# print(is_short_chain_fatty_acid("CCCC(O)=O"))  # Represents butyric acid