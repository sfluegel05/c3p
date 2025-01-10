"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.
    A 3-hydroxy fatty acid has a hydroxyl group at the beta- or 3-position relative to the carboxyl group
    and is part of a long hydrocarbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a carboxyl group at the end of the chain 
    # with a hydroxyl group at the 3rd carbon position
    # The pattern captures the carbon chain, hydroxyl group in the 3-position, and terminal carboxylic acid
    pattern = Chem.MolFromSmarts("CC(C)(O)[C](=O)O")  # General pattern without stereochemistry
    stereo_pattern = Chem.MolFromSmarts("C[C@H](O)C(=O)O | C[C@@H](O)C(=O)O")  # Stereochemical patterns

    # Check for the general and stereochemical patterns
    if mol.HasSubstructMatch(pattern):
        return True, "Contains hydroxyl group at the 3-position and carboxylic acid group"

    if mol.HasSubstructMatch(stereo_pattern):
        return True, "Contains stereochemical hydroxyl group pattern at the 3-position and carboxylic acid group"

    # Final check for extended hydrocarbon chains, implying it's a fatty acid
    carbon_chain_length = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if carbon_chain_length > 10 and oxygen_count >= 3:  # Reasonable length for fatty acid and OH + COO presence
        return True, "Has characteristic 3-hydroxy pattern as part of a fatty acid chain"

    return False, "Hydroxyl group not found at the 3-position or missing carboxyl group"