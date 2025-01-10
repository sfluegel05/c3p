"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
from rdkit import Chem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    An alkanesulfonate oxoanion contains a sulfonate group attached to an aliphatic carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanesulfonate oxoanion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define sulfonate group pattern, potentially considering different variants
    sulfonate_pattern = Chem.MolFromSmarts("[S](=O)(=O)[O-]")  
    if not mol.HasSubstructMatch(sulfonate_pattern):
        return False, "Missing sulfonate group S(=O)(=O)[O-]"

    # Initialize flag to detect aliphatic attachment to sulfonate
    aliphatic_connected = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'S':
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C' and not neighbor.GetIsAromatic():
                    # Check for linear aliphatic chain length, seeking chains of at least three carbons
                    chain_length = 0
                    current_atom = neighbor
                    visited = set()
                    stack = [current_atom]
                    while stack:
                        current = stack.pop()
                        if current.GetIdx() in visited:
                            continue
                        visited.add(current.GetIdx())

                        if current.GetSymbol() == 'C' and not current.GetIsAromatic():
                            chain_length += 1
                            for next_atom in current.GetNeighbors():
                                if next_atom.GetIdx() not in visited and next_atom.GetSymbol() == 'C':
                                    stack.append(next_atom)

                    if chain_length >= 2:
                        aliphatic_connected = True
                        break
    
    if aliphatic_connected:
        return True, "Sulfonate group attached to an aliphatic structure"
        
    return False, "Sulfonate group not connected to an adequate aliphatic structure"