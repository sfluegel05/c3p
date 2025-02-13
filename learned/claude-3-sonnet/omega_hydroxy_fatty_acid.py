"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
"""
Classifies: omega-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an omega-hydroxy fatty acid based on its SMILES string.
    An omega-hydroxy fatty acid has a carboxyl group at position 1 and a hydroxyl at the
    opposite end (omega position) of a primarily linear carbon chain.
    
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

    # Get the indices of key atoms
    carboxylic_carbon = mol.GetSubstructMatch(carboxylic_pattern)[0]
    terminal_ch2 = mol.GetSubstructMatch(terminal_oh_pattern)[0]
    
    # Check for large rings (allowing epoxides)
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) > 3:  # Allow 3-membered rings (epoxides)
            return False, "Contains non-epoxide ring structure"
    
    # Get shortest path between COOH carbon and terminal CH2
    path = Chem.GetShortestPath(mol, carboxylic_carbon, terminal_ch2)
    if not path:
        return False, "No valid path between COOH and terminal OH"
    
    chain_length = len(path)
    if chain_length < 3:
        return False, "Chain too short for fatty acid"
    
    # Count atoms to verify composition
    total_atoms = mol.GetNumAtoms()
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Check carbon content (more permissive than before)
    if carbon_count < total_atoms * 0.4:
        return False, "Not primarily carbon-based"
    
    # Check oxygen content (must have at least 3: COOH + OH)
    if oxygen_count < 3:
        return False, "Insufficient oxygen atoms for omega-hydroxy fatty acid"
    
    # Check for excessive branching by comparing path length to total carbons
    # Allow more branching for longer chains
    max_allowed_extra_carbons = max(chain_length // 3, 2)  # More permissive
    if carbon_count > chain_length + max_allowed_extra_carbons:
        return False, "Too many branches from main chain"
    
    # Verify the molecule consists mainly of C, H, O
    allowed_atoms = {1, 6, 8}  # H, C, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, "Contains elements other than C, H, O"

    return True, "Contains terminal hydroxyl and carboxylic acid groups at opposite ends of primarily linear carbon chain"