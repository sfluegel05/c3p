"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
"""
Classifies: omega-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.GraphDescriptors import ComplexityPickler
from collections import deque

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
    terminal_ch2_oh = mol.GetSubstructMatch(terminal_oh_pattern)[0]
    
    # Check for rings (allowing small epoxides)
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) > 3:  # Allow 3-membered rings (epoxides)
            return False, "Contains non-epoxide ring structure"
            
    # Find path between COOH and CH2OH
    paths = Chem.FindAllPathsOfLengthN(mol, 
                                     min(mol.GetNumAtoms()-1, 30),  # Limit path length
                                     carboxylic_carbon,
                                     terminal_ch2_oh)
    
    if not paths:
        return False, "No valid path between COOH and terminal OH"
    
    # Get the shortest path
    main_chain = min(paths, key=len)
    chain_length = len(main_chain)
    
    if chain_length < 3:
        return False, "Chain too short for fatty acid"
        
    # Check that the path is primarily carbons
    non_carbon_count = sum(1 for idx in main_chain 
                         if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6)
    if non_carbon_count > chain_length * 0.3:  # Allow up to 30% non-carbon
        return False, "Main chain not primarily carbon"
        
    # Count branches along main chain
    branch_count = 0
    for idx in main_chain:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 6:  # Carbon atom
            # Count carbon neighbors not in main chain
            c_neighbors = sum(1 for n in atom.GetNeighbors() 
                            if n.GetAtomicNum() == 6 and n.GetIdx() not in main_chain)
            if c_neighbors > 1:  # Allow one branch
                return False, "Too many branches from main chain"
            branch_count += c_neighbors
            
    if branch_count > chain_length/3:
        return False, "Too many total branches for fatty acid"
        
    # Check that terminal groups are actually at the ends of the main chain
    if carboxylic_carbon != main_chain[0] and carboxylic_carbon != main_chain[-1]:
        return False, "Carboxylic group not at end of main chain"
    if terminal_ch2_oh != main_chain[0] and terminal_ch2_oh != main_chain[-1]:
        return False, "Terminal hydroxyl not at end of main chain"
        
    # Get all atoms and check overall composition
    total_atoms = len(mol.GetAtoms())
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # For very short chains (<=4 carbons), use different thresholds
    if chain_length <= 4:
        if carbon_count < total_atoms * 0.4:  # More permissive for short chains
            return False, "Not primarily carbon-based (short chain)"
    else:
        if carbon_count < total_atoms * 0.6:  # Stricter for longer chains
            return False, "Not primarily carbon-based"
    
    if oxygen_count < 3:
        return False, "Insufficient oxygen atoms for omega-hydroxy fatty acid"
    if oxygen_count > carbon_count/2:
        return False, "Too many oxygen atoms for fatty acid structure"

    return True, "Contains terminal hydroxyl and carboxylic acid groups at opposite ends of linear carbon chain"