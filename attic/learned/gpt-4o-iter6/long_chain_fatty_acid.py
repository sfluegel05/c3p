"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acid based on its SMILES string.
    A long-chain fatty acid is defined as a fatty acid with a carbon chain length from C13 to C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylic acid group (C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Get all carbon atoms
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    
    # Calculate chain length (excluding carboxylic acid carbon and neighboring oxygen)
    # Match MCC => carboxylic acid carbon, which is terminal and follow an oxygen, remove it 
    for atom in carbon_atoms:
        if mol.GetAtomWithIdx(atom.GetIdx()).GetDegree() == 3:  # The degree attribute refers to the number of directly bonded atoms
            carbon_atoms.remove(atom)
            break
    
    chain_length = len(carbon_atoms)
    
    if chain_length >= 13 and chain_length <= 22:
        return True, f"Valid long chain fatty acid with {chain_length} carbons"
    else:
        return False, f"Chain length {chain_length} out of range for long-chain fatty acids (C13-C22)"

# Example usage
smiles_example = "C(CCCCCCCCCCCCCCC)C(O)=O"
result, reason = is_long_chain_fatty_acid(smiles_example)
print(result, reason)