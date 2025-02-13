"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
"""
Classifies: omega-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an omega-hydroxy fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an omega-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
    
    # Count carboxylic acid groups - should only have one
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_pattern)
    if len(carboxylic_matches) > 1:
        return False, "Multiple carboxylic acid groups found"
    
    # Look for terminal hydroxyl group (CH2-OH)
    terminal_oh_pattern = Chem.MolFromSmarts("[CH2][OX2H1]")
    if not mol.HasSubstructMatch(terminal_oh_pattern):
        return False, "No terminal hydroxyl group found"
    
    # Get all carbons and oxygens
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Check minimum chain length (at least 3 carbons)
    if carbon_count < 3:
        return False, "Chain too short for fatty acid"
    
    # Verify it's primarily a carbon chain
    # Most atoms should be carbons for a fatty acid
    total_atoms = len(mol.GetAtoms())
    if carbon_count < total_atoms * 0.6:  # At least 60% should be carbon
        return False, "Not primarily a carbon chain"
        
    # Check for reasonable number of oxygens
    # Must have at least 3 (COOH + OH), but not too many relative to chain length
    if oxygen_count < 3:
        return False, "Insufficient oxygen atoms for omega-hydroxy fatty acid"
    if oxygen_count > carbon_count/2:
        return False, "Too many oxygen atoms for fatty acid structure"
    
    # Verify basic connectivity
    # Get the carboxylic carbon atom
    carboxylic_match = mol.GetSubstructMatch(carboxylic_pattern)
    if not carboxylic_match:
        return False, "Cannot find carboxylic group connectivity"
    
    carboxylic_carbon = mol.GetAtomWithIdx(carboxylic_match[0])
    
    # Check if molecule is too branched to be a fatty acid
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            if len([n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]) > 2:
                # More than two carbon neighbors means significant branching
                return False, "Structure too branched for fatty acid"

    return True, "Contains terminal hydroxyl and carboxylic acid group in linear chain structure"