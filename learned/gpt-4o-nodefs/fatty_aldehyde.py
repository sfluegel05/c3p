"""
Classifies: CHEBI:35746 fatty aldehyde
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    Fatty aldehydes typically have long carbon chains and terminal aldehyde groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty aldehyde, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for terminal aldehyde group: formaldehyde structure `C=O` connected to hydrogen
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[CH3]")
    if not mol.HasSubstructMatch(aldehyde_pattern):
        # Check an alternate pattern ensuring C=O but being flexible
        aldehyde_pattern_alt = Chem.MolFromSmarts("[CX3H1](=O)[CH2]")
        if not mol.HasSubstructMatch(aldehyde_pattern_alt):
            return False, "No terminal aldehyde group found"

    # Determine the main carbon chain length
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Assume fatty aldehydes to have a minimum chain length, but sometimes higher than 5
    min_chain_length = 6  # Adjust higher to accommodate typical examples
    
    if carbon_count < min_chain_length:
        return False, f"Carbon chain too short for fatty aldehyde, found {carbon_count} carbons"
    
    # Unsaturation presence is typical but not a hard rule
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    has_double_bonds = mol.HasSubstructMatch(double_bond_pattern)
    
    return True, f"Contains terminal aldehyde group and a carbon chain with {carbon_count} carbons. Unsaturated: {'Yes' if has_double_bonds else 'No'}"