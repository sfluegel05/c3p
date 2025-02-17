"""
Classifies: CHEBI:18000 aralkylamine
"""
"""
Classifies: aralkylamine
Definition: An alkylamine in which one or more alkyl substituents bears an aromatic group.
Operationally, we require that at least one non‐aromatic amine nitrogen (i.e. not in an aromatic ring)
has a branch (a substituent) that contains an aromatic carbon atom (a marker for an aromatic ring)
within 1 to 3 bonds from that nitrogen.
"""

from rdkit import Chem
from collections import deque

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    
    The approach is:
      1. Parse the SMILES.
      2. Identify non‐aromatic nitrogen atoms (atomic number 7 which are not flagged as aromatic).
      3. For each such nitrogen, traverse its substituent branches (using a breadth‐first search,
         up to a bond distance of 3) and check if an aromatic carbon (atomic number 6 flagged as aromatic)
         is encountered.
      4. If found, return True with the bond-distance at which the aromatic carbon was encountered.
    
    Args:
      smiles (str): SMILES string for the molecule.
    
    Returns:
      bool: True if molecule is classified as an aralkylamine, False otherwise.
      str: Explanation (reason) for classification.
      
    If the SMILES cannot be parsed, return (False, "Invalid SMILES string").
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all non-aromatic nitrogen atoms
    non_aromatic_nitrogens = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7 and not atom.GetIsAromatic()]
    
    if not non_aromatic_nitrogens:
        return False, "No non‐aromatic amine nitrogen found in the molecule."
    
    max_depth = 3  # Maximum bond separation allowed
    
    # For each non aromatic nitrogen, search its substituent branch
    for n_atom in non_aromatic_nitrogens:
        n_idx = n_atom.GetIdx()
        
        # We'll use a breadth-first search from this nitrogen.
        # Each entry in the queue is (current_atom_index, depth) where depth is the number of bonds from the nitrogen.
        visited = set([n_idx])
        queue = deque()
        
        # Start from immediate neighbors of the nitrogen (depth=1)
        for neighbor in n_atom.GetNeighbors():
            neigh_idx = neighbor.GetIdx()
            if neigh_idx not in visited:
                queue.append((neigh_idx, 1))
                visited.add(neigh_idx)
                
        while queue:
            current_idx, depth = queue.popleft()
            current_atom = mol.GetAtomWithIdx(current_idx)
            # Check: if current atom is an aromatic carbon then we have an aromatic substituent
            if current_atom.GetAtomicNum() == 6 and current_atom.GetIsAromatic():
                return True, (f"Found a non‐aromatic amine nitrogen (atom {n_idx}) with an aromatic substituent "
                              f"(atom {current_idx}) at a bond distance of {depth}. Molecule classified as aralkylamine.")
            
            # Only continue searching if we haven't reached maximum depth.
            if depth < max_depth:
                for neighbor in current_atom.GetNeighbors():
                    neigh_idx = neighbor.GetIdx()
                    if neigh_idx not in visited:
                        queue.append((neigh_idx, depth + 1))
                        visited.add(neigh_idx)
                        
    return False, ("No aralkylamine substructure found: no non‐aromatic amine is connected "
                   "to an aromatic carbon (marker for aromatic group) within 3 bonds in its substituent branch.")

# Optional test block; remove or comment out when deploying as a module.
if __name__ == "__main__":
    # Some test cases based on our working definition.
    test_smiles = [
        ("NCc1ccccc1", "benzylamine"),           # Expected: True (distance 2)
        ("NCCc1ccccc1", "phenethylamine"),         # Expected: True (distance 3)
        ("c1ccc(N)cc1", "aniline (should not classify, N aromatic)"),   # Expected: False
        ("OCCNC1=CC=CC=C1", "2-Anilinoethanol")     # Expected: True
    ]
    
    for sm, name in test_smiles:
        result, reason = is_aralkylamine(sm)
        print(f"SMILES: {sm}\nMolecule: {name}\nClassification: {result}\nReason: {reason}\n")