"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdchem

def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17beta-hydroxy steroid based on its SMILES string.
    A 17beta-hydroxy steroid has the basic steroid core and a beta-configured hydroxyl group at position 17.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for the steroid core
    steroid_core_pattern = Chem.MolFromSmarts("[C]1[C][C][C]2[C][C]([C])([C][C]C)[C][C]2[C]1")
    
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Molecule does not have the basic steroid core structure."
    
    # SMARTS pattern for C17
    c17_pattern = Chem.MolFromSmarts("[C]1~[C]~[C]~[#6,#7,#8,#9]")
    c17_matches = mol.GetSubstructMatches(c17_pattern)

    if not c17_matches:
       return False, "C17 not found"
    
    c17_idx = None
    for match in c17_matches:
      c17_idx = match[0]
      c17_atom = mol.GetAtomWithIdx(c17_idx)
      
      # Verify that this C is actually at a steroid core
      is_c17_in_core = False
      for core_match in mol.GetSubstructMatches(steroid_core_pattern):
        for core_atom_idx in core_match:
            core_atom = mol.GetAtomWithIdx(core_atom_idx)
            if c17_atom.GetIdx() == core_atom.GetIdx():
               is_c17_in_core = True
               break
        if is_c17_in_core:
           break
      if is_c17_in_core:
          break
    if not is_c17_in_core:
          return False, "C17 not in a steroid core"

    # Check for hydroxyl group at C17
    has_oh = False
    for neighbor in c17_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 8: # Oxygen atom
          has_oh = True
          oh_idx = neighbor.GetIdx()
          break
    if not has_oh:
       return False, "No hydroxyl group found at C17"

    # Check the chirality of C17.
    chiral_tag = c17_atom.GetChiralTag()
    if chiral_tag == rdchem.ChiralType.CHI_TETRAHEDRAL_CW or chiral_tag == rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
       
       # Now verify that it's a beta hydroxyl group.
       # Get the neighboring atoms of C17 other than O
       neighbor_atoms = []
       for nbr in c17_atom.GetNeighbors():
          if nbr.GetIdx() != oh_idx:
              neighbor_atoms.append(nbr)
       
       # If 3 neighbors
       if len(neighbor_atoms) == 3:
         # Find a carbon substituent
         c_idx_1 = None
         c_idx_2 = None

         for nbr in neighbor_atoms:
            if nbr.GetAtomicNum() == 6:
                if c_idx_1 is None:
                  c_idx_1 = nbr.GetIdx()
                else:
                  c_idx_2 = nbr.GetIdx()
         if c_idx_1 is None or c_idx_2 is None:
              return False, "No two carbon neighbours around C17"
         c1_atom = mol.GetAtomWithIdx(c_idx_1)
         c2_atom = mol.GetAtomWithIdx(c_idx_2)

         # Get bond order with one of the two carbons.
         bond_order_c17_c1 = mol.GetBondBetweenAtoms(c17_idx, c_idx_1).GetBondType()
         bond_order_c17_c2 = mol.GetBondBetweenAtoms(c17_idx, c_idx_2).GetBondType()

         # Check the chirality with the 2 carbon atoms.
         if bond_order_c17_c1 == Chem.BondType.SINGLE and bond_order_c17_c2 == Chem.BondType.SINGLE:
           # Use RDKit to check the chirality of the C17. We look at the positions of the neighbors to check that the OH is at the correct position (beta)
           if chiral_tag == rdchem.ChiralType.CHI_TETRAHEDRAL_CW:
                # Clockwise (beta, going into the plane of the paper)
                
              conf = mol.GetConformer()
              c17_pos = conf.GetAtomPosition(c17_idx)
              c1_pos = conf.GetAtomPosition(c_idx_1)
              c2_pos = conf.GetAtomPosition(c_idx_2)
              oh_pos = conf.GetAtomPosition(oh_idx)

              v1 = (c1_pos - c17_pos).normalized()
              v2 = (c2_pos - c17_pos).normalized()
              v3 = (oh_pos - c17_pos).normalized()

              normal = v1.CrossProduct(v2)

              if normal.DotProduct(v3) < 0:
                 return False, "Chirality does not match 17-beta"
           elif chiral_tag == rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
                # Counterclockwise
              conf = mol.GetConformer()
              c17_pos = conf.GetAtomPosition(c17_idx)
              c1_pos = conf.GetAtomPosition(c_idx_1)
              c2_pos = conf.GetAtomPosition(c_idx_2)
              oh_pos = conf.GetAtomPosition(oh_idx)

              v1 = (c1_pos - c17_pos).normalized()
              v2 = (c2_pos - c17_pos).normalized()
              v3 = (oh_pos - c17_pos).normalized()

              normal = v1.CrossProduct(v2)
              if normal.DotProduct(v3) > 0:
                  return False, "Chirality does not match 17-beta"
         else:
             return False, "C17 has non-single carbon neighbors"
    elif chiral_tag == rdchem.ChiralType.CHI_UNSPECIFIED or chiral_tag == rdchem.ChiralType.CHI_NONE:
       # If chirality not specified then it must be that the only substituent is H. If not, it's not beta.
       count_h = 0
       for neighbor in c17_atom.GetNeighbors():
          if neighbor.GetAtomicNum() == 1:
             count_h += 1
       if count_h != 1:
        return False, "Chirality must be specified for non-H 17 substituent"
    else:
       return False, "Unknown chiral tag."
       
    return True, "Molecule contains a steroid core with a beta-hydroxy group at position 17."