"""
Classifies: CHEBI:50753 isoflavonoid
"""
"""
Classifies: isoflavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    An isoflavonoid is defined as any 1-benzopyran with an aryl substituent at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the chromone core with a variable substituent at position 3
    # The SMARTS pattern includes atom mapping to identify position 3
    chromone_smarts = '[O]=C1C=CC([#6:1])=CC=C1'  # Atom map 1 corresponds to position 3
    pattern = Chem.MolFromSmarts(chromone_smarts)
    if pattern is None:
        return False, "Invalid SMARTS pattern for chromone core"
    
    # Build mapping from query atom indices to atom map numbers
    query_map = {}  # atom map number to query atom index
    for atom in pattern.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if map_num > 0:
            query_map[map_num] = atom.GetIdx()
    
    # Find substructure matches of the chromone core
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "No chromone core found"
    
    for match in matches:
        # Get the atom index corresponding to position 3
        atom_idx = match[query_map[1]]
        atom = mol.GetAtomWithIdx(atom_idx)
        
        # Get the substituents attached to position 3 that are not part of the chromone core
        pattern_atom_indices = set(match)
        neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() not in pattern_atom_indices]
        if not neighbors:
            continue  # No substituents attached to position 3
        
        # Check if the substituent contains an aryl group
        aryl_found = False
        for nbr in neighbors:
            # Perform a breadth-first search to find an aromatic ring
            visited = set()
            stack = [nbr]
            while stack:
                current_atom = stack.pop()
                current_idx = current_atom.GetIdx()
                if current_idx in visited or current_idx in pattern_atom_indices:
                    continue
                visited.add(current_idx)
                if current_atom.GetIsAromatic() and current_atom.IsInRing():
                    aryl_found = True
                    break
                else:
                    stack.extend([n for n in current_atom.GetNeighbors() if n.GetIdx() not in visited])
            if aryl_found:
                break  # Aryl group found
        if aryl_found:
            return True, "Contains chromone core with aryl substituent at position 3"
        else:
            continue  # Check next match if available
    
    return False, "Chromone core found, but no aryl substituent at position 3"

# Example usage:
# smiles = 'O=C1C=CC2=CC=CC=C2O1C3=CC=CC=C3'  # Isoflavone
# result, reason = is_isoflavonoid(smiles)
# print(result, reason)