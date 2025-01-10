"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
from rdkit import Chem

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an omega-hydroxy fatty acid based on its SMILES string.
    An omega-hydroxy fatty acid is characterized by a carboxyl group at position 1 and a hydroxyl group at the omega (last) position.

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

    # Get list of atoms
    atoms = mol.GetAtoms()
    
    # Identify carboxyl group ('C(=O)O') at the beginning
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No terminal carboxyl group found"
    
    # Track beginning of the chain
    carboxyl_end = carboxyl_matches[0][-1]  # Last atom of the pattern is the O in 'C(=O)O'

    # Identify hydroxyl group (OH) at the end opposite the carboxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("O")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    # Check if any hydroxyl is really terminal and is opposite to carboxyl
    terminal_hydroxyl = False
    for match in hydroxyl_matches:
        # Check if the hydroxyl group is at the chain's end
        atom_idx = match[0]
        atom = mol.GetAtomWithIdx(atom_idx)
        neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
        
        # For a terminal hydroxyl at the omega position, it should have no carbon neighbor further than one step from non-hydroxyl atoms 
        if len(neighbors) == 1 and atom.GetDegree() == 1:
            chain_length_through_oxygen = mol.GetDistanceMatrix()[carboxyl_end][atom_idx]
            if chain_length_through_oxygen >= 5:  # Minimum chain length between carboxyl and hydroxyl
                terminal_hydroxyl = True
                break
    
    if not terminal_hydroxyl:
        return False, "No terminal omega-hydroxyl group found"

    # Ensure all carbon atoms are in one primary chain for a straight configuration
    carbon_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
    if carbon_count < 6:
        return False, f"Chain length too short for fatty acid (found {carbon_count} carbons)"
    if carboxyl_end >= carbon_count:
        return False, "Carboxyl at non-primary chain position"

    return True, "Matches omega-hydroxy fatty acid structure"