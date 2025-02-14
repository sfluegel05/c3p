"""
Classifies: CHEBI:15889 sterol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.
    A sterol is a steroid with a hydroxyl group at the 3-beta position.
     It has a tetracyclic ring system closely related to cholestan-3-ol (additional carbon atoms may be present in the side chain).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Steroid core SMARTS pattern (tetracyclic ring system with 6-6-6-5 fusion)
    # This pattern is more generic, focusing on the ring structure.
    steroid_core_pattern = Chem.MolFromSmarts("[C]1[C]2[C]3[C]([C]1)[C][C]4[C]([C]3[C]2)[C][C]5[C]([C]4)[C][C]5")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "No steroid core found"


    core_match = mol.GetSubstructMatch(steroid_core_pattern)
    if not core_match:
        return False, "No match of steroid_core_pattern"


    # Identify C3 atom - this may not be carbon 2 every time.  We now find the atoms that form the 4 rings, and check for the hydroxyl group on an atom in one of the outer rings.
    c3_candidates = []
    
    #We will have to identify our c3 by checking for a neighbor carbon that has two neighbors.
    ring_atoms_idx = core_match

    for atom_idx in ring_atoms_idx:
        c_atom = mol.GetAtomWithIdx(atom_idx)
        
        if c_atom.GetSymbol() == 'C':
          neighbors = c_atom.GetNeighbors()
          if len(neighbors) == 3 and  all(n.GetSymbol() == "C" for n in neighbors):
            c3_candidates.append(c_atom)

    # Check for 3-beta hydroxyl group. Must be attached to a candidate of the 3 position and have one H.
    has_3beta_oh = False
    for c3 in c3_candidates:
      for neighbor in c3.GetNeighbors():
          if neighbor.GetSymbol() == 'O' and neighbor.GetTotalNumHs() == 1:
            has_3beta_oh = True
            break
      if has_3beta_oh:
        break
    
    if not has_3beta_oh:
      return False, "No hydroxyl group at the 3-beta position"
      
    
    # Identify a potential C17 atom.  We find an atom that has three other neighbors and two of those neighbors have 2 or less neighbors.
    c17_candidates = []
    for atom_idx in ring_atoms_idx:
       c_atom = mol.GetAtomWithIdx(atom_idx)
       if c_atom.GetSymbol() == 'C':
          neighbors = c_atom.GetNeighbors()
          if len(neighbors) == 3:
             num_neighbors_two = sum([1 if len(n.GetNeighbors()) <=2 else 0 for n in neighbors])
             if num_neighbors_two >= 2:
              c17_candidates.append(c_atom)

    # Check for side chain at a potential C17 (at least one carbon atom attached)
    has_sidechain = False
    for c17 in c17_candidates:
       for neighbor in c17.GetNeighbors():
           if neighbor.GetSymbol() == 'C':
               has_sidechain = True
               break
       if has_sidechain:
          break
    if not has_sidechain:
        return False, "No side chain at C17"


    return True, "Contains a steroid core, a hydroxyl group at the 3-beta position, and a side chain at C17."