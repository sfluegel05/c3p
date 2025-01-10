"""
Classifies: CHEBI:35746 fatty aldehyde
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    Fatty aldehydes typically have long carbon chains and terminal aldehyde groups (R-CHO).

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
    
    # Look for terminal aldehyde group: formaldehyde-like structure `H-C=O`
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[CX4H3]")
    if not mol.HasSubstructMatch(aldehyde_pattern):
        return False, "No terminal aldehyde group found"

    # Determine the main carbon chain length
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Set a minimum chain length typical for fatty aldehydes
    min_chain_length = 8  # Adjust based on known fatty aldehyde examples
    
    if carbon_count < min_chain_length:
        return False, f"Carbon chain too short for fatty aldehyde, found {carbon_count} carbons"
    
    # Unsaturation presence: checked but not necessary for the definition
    has_double_bonds = any(bond.GetBondType() == Chem.rdchem.BondType.DOUBLE for bond in mol.GetBonds())
    
    return True, f"Contains terminal aldehyde group and a carbon chain with {carbon_count} carbons. Unsaturated: {'Yes' if has_double_bonds else 'No'}"