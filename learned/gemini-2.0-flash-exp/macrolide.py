"""
Classifies: CHEBI:25106 macrolide
"""
from rdkit import Chem
from collections import deque

def is_macrolide(smiles: str):
    """
    Determines if a molecule is a macrolide based on its SMILES string.
    Macrolides are characterized by a macrocyclic lactone ring (12+ atoms) derived from polyketides.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macrolide, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find lactone group
    lactone_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2]")
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)

    if not lactone_matches:
          return False, "No lactone group found"
    
    for match in lactone_matches: # check all lactone groups until we find a macrolide

        lactone_carbon_index = match[0]

        # Breadth-first search to find the ring
        queue = deque([(lactone_carbon_index, [lactone_carbon_index])])  # (atom_index, path)
        visited_atoms = set() #to track the ring size

        while queue:
            current_atom_index, path = queue.popleft()
            
            if current_atom_index in visited_atoms:
                continue
            visited_atoms.add(current_atom_index)
            
            if len(visited_atoms) > 1:
                if current_atom_index == lactone_carbon_index:
                    if len(visited_atoms) >= 12:
                         return True, f"Macrocyclic lactone ring (size: {len(visited_atoms)}) found."
                    else:
                       break # continue the loop, to check the other lactone
            
            
            atom = mol.GetAtomWithIdx(current_atom_index)
            for neighbor in atom.GetNeighbors():
                 neighbor_index = neighbor.GetIdx()
                 if neighbor_index not in path:
                     queue.append((neighbor_index, path+[neighbor_index]))
                    
    return False, "No macrocyclic lactone ring (12+ atoms) found."