"""
Classifies: CHEBI:35746 fatty aldehyde
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty aldehyde, False otherwise
        str: Reason for classification
    """
    
    # Try to parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for aldehyde group C=O with terminal H
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")  # Terminal aldehyde carbon
    if not mol.HasSubstructMatch(aldehyde_pattern):
        return False, "No terminal aldehyde group found"
    
    # Count carbons in the molecule
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Check for presence of at least a moderate carbon chain
    # Adjust threshold based on general knowledge about fatty aldehyde chain length
    min_chain_length = 5  # Typical chain length might be longer, adjust empirically
    
    if carbon_count < min_chain_length:
        return False, f"Carbon chain too short for fatty aldehyde, found {carbon_count} carbons"
    
    # Check for presence of unsaturation (double bonds)
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(double_bond_pattern):
        return False, "No carbon-carbon double bonds found (not required, but typical)"
    
    return True, f"Contains terminal aldehyde group and a carbon chain with {carbon_count} carbons"

# Example usage
smiles_examples = [
    "O=C(/C=C/C(=O)C=O)C", "O=C/C=C/C=C\CCC", "CCCCCCCC=O", "O=CCC\C=C\CCCC"
]
for smiles in smiles_examples:
    is_fatty, reason = is_fatty_aldehyde(smiles)
    print(f"SMILES: {smiles}, Is Fatty Aldehyde: {is_fatty}, Reason: {reason}")