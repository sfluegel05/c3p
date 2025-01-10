"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
"""
Classifies: CHEBI:XXXXX omega-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an omega-hydroxy fatty acid based on its SMILES string.
    An omega-hydroxy fatty acid is a straight-chain fatty acid with a carboxyl group at position 1
    and a hydroxyl group at the terminal (omega) position.

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

    # Check for carboxyl group (-COOH) at position 1
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No carboxyl group found"

    # Find the longest carbon chain
    # Use a more flexible pattern to handle double bonds and rings
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if not chain_matches:
        return False, "No long carbon chain found"

    # Find the terminal carbon of the longest chain
    terminal_carbon = None
    for match in chain_matches:
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetDegree() == 1 or (atom.GetDegree() == 2 and atom.GetAtomicNum() == 6):  # Terminal carbon can have degree 1 or 2 (if part of a double bond)
                terminal_carbon = atom_idx
                break
        if terminal_carbon is not None:
            break

    if terminal_carbon is None:
        return False, "No terminal carbon found in the longest chain"

    # Check if the terminal carbon has a hydroxyl group
    terminal_atom = mol.GetAtomWithIdx(terminal_carbon)
    has_hydroxyl = False
    for neighbor in terminal_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() == 1:  # Oxygen with one hydrogen (OH)
            has_hydroxyl = True
            break

    if not has_hydroxyl:
        return False, "No hydroxyl group found at the terminal position"

    # Check if the molecule is a straight-chain fatty acid (no branching)
    # Count the number of carbons in the longest chain
    chain_length = len(chain_matches[0])
    if chain_length < 6:
        return False, "Chain too short to be a fatty acid"

    # Ensure no branching in the chain
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetDegree() > 2:  # Carbon with more than 2 non-hydrogen connections (branching)
            return False, "Branched chain detected"

    return True, "Contains a carboxyl group at position 1 and a hydroxyl group at the terminal position of a straight-chain fatty acid"