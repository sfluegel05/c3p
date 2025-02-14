"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
"""
Classifies: omega-hydroxy fatty acids
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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
    carboxylic_pattern = Chem.MolFromSmarts("[OX2H][CX3](=[OX1])")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
    
    # Find all hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    
    if len(hydroxyl_matches) < 2:  # Need at least 2 (one from COOH, one terminal)
        return False, "Insufficient hydroxyl groups"
    
    # Get atoms for analysis
    atoms = mol.GetAtoms()
    
    # Find terminal carbons (connected to only one other carbon)
    terminal_carbons = []
    for atom in atoms:
        if atom.GetAtomicNum() == 6:  # Carbon atom
            carbon_neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
            if len(carbon_neighbors) == 1:
                terminal_carbons.append(atom.GetIdx())
    
    # Check if any terminal carbon has a hydroxyl group (excluding the carboxylic acid carbon)
    terminal_oh_found = False
    carboxyl_carbon = mol.GetSubstructMatch(carboxylic_pattern)[1]  # Index of C in COOH
    
    for tc in terminal_carbons:
        if tc != carboxyl_carbon:  # Skip the carboxylic acid carbon
            atom = mol.GetAtomWithIdx(tc)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() == 1:
                    terminal_oh_found = True
                    break
    
    if not terminal_oh_found:
        return False, "No terminal hydroxyl group found"
    
    # Check for reasonable chain length (at least 3 carbons)
    carbon_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
    if carbon_count < 3:
        return False, "Carbon chain too short"
    
    # Verify it's primarily a chain structure (can have branches but should be mostly linear)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 1:  # Allow up to 1 ring (for epoxy groups)
        return False, "Too many ring structures"
    
    # Check if the structure is primarily linear (most carbons should have 1-2 carbon neighbors)
    linear_carbons = 0
    for atom in atoms:
        if atom.GetAtomicNum() == 6:
            carbon_neighbors = len([n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6])
            if carbon_neighbors <= 2:
                linear_carbons += 1
    
    if linear_carbons < carbon_count * 0.7:  # At least 70% should be part of linear chain
        return False, "Structure not sufficiently linear"
        
    return True, "Contains terminal hydroxyl and carboxylic acid groups in a primarily linear structure"