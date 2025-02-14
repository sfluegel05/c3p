"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy-Delta(5)-steroid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3beta-hydroxy-Delta(5)-steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for the steroid core (4 fused rings, 3 six membered and 1 five membered)
    steroid_core_pattern = Chem.MolFromSmarts("[r6]1[r6][r6][r5]1")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Steroid core not found"

    # Define SMARTS pattern for the 3-beta-hydroxyl group, without strict chirality.
    beta_hydroxy_pattern = Chem.MolFromSmarts("[CX4][OH1]")
    hydroxy_matches = mol.GetSubstructMatches(beta_hydroxy_pattern)
    
    has_3beta_hydroxy = False
    for match in hydroxy_matches:
        hydroxy_carbon = match[0]
        #check if the carbon with the hydroxyl is part of the steroid core
        is_on_steroid = False
        for core_match in mol.GetSubstructMatches(steroid_core_pattern):
            for core_atom_idx in core_match:
                if core_atom_idx == hydroxy_carbon:
                    is_on_steroid = True
                    break
            if is_on_steroid:
                break
        if is_on_steroid:
            #check if it is the "3-beta" by checking that it is connected to carbon of the core
            # In steroids C3 is bonded to 2 other carbons of the core, so we verify that the hydroxy carbon is bonded to two core carbons
            core_carbons_count = 0
            for neighbor in mol.GetAtomWithIdx(hydroxy_carbon).GetNeighbors():
               for core_match in mol.GetSubstructMatches(steroid_core_pattern):
                 if neighbor.GetIdx() in core_match:
                      core_carbons_count += 1
            if core_carbons_count == 2:
              has_3beta_hydroxy = True
              break

    if not has_3beta_hydroxy:
       return False, "No 3-beta hydroxyl group found on the steroid core"
    
    # Define SMARTS pattern for the double bond between C5 and C6
    double_bond_pattern = Chem.MolFromSmarts("[C]=[C]")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    
    has_delta5_bond = False
    for match in double_bond_matches:
        carbon1_idx = match[0]
        carbon2_idx = match[1]

        is_c5_c6 = False
        # Check if both carbons are part of the steroid core
        is_c1_in_core = False
        is_c2_in_core = False
        for core_match in mol.GetSubstructMatches(steroid_core_pattern):
            if carbon1_idx in core_match:
                is_c1_in_core = True
            if carbon2_idx in core_match:
               is_c2_in_core = True
        if not (is_c1_in_core and is_c2_in_core):
             continue
       
        # Check if the bond is between C5 and C6
        # We will count the number of carbons separating each of the double bond carbons to a specific corner carbon (C10).
        # C5 is typically 2 bonds away from C10 and C6 is 1 bond away, so we will verify that one carbon is 1 away from a corner and the other is 2 away.
        # There are multiple "corners", so we pick the carbons that are linked to 3 other core carbons as "corners" and do a search for both of them.
        corner_atoms = []
        for core_match in mol.GetSubstructMatches(steroid_core_pattern):
            for core_atom_idx in core_match:
                core_atom = mol.GetAtomWithIdx(core_atom_idx)
                core_carbons_count = 0
                for neighbor in core_atom.GetNeighbors():
                  if neighbor.GetIdx() in core_match:
                     core_carbons_count += 1
                if core_carbons_count == 3:
                    corner_atoms.append(core_atom_idx)
        if len(corner_atoms) < 2 :
              continue
        for corner1 in corner_atoms:
          for corner2 in corner_atoms:
            if corner1 == corner2:
                continue

            # perform bond count
            bond_count_1_corner1 = len(Chem.GetShortestPath(mol,carbon1_idx,corner1)) - 1
            bond_count_2_corner1 = len(Chem.GetShortestPath(mol,carbon2_idx,corner1)) - 1
            bond_count_1_corner2 = len(Chem.GetShortestPath(mol,carbon1_idx,corner2)) - 1
            bond_count_2_corner2 = len(Chem.GetShortestPath(mol,carbon2_idx,corner2)) - 1

            if (bond_count_1_corner1 == 2 and bond_count_2_corner1 == 1) or \
               (bond_count_1_corner1 == 1 and bond_count_2_corner1 == 2) or \
               (bond_count_1_corner2 == 2 and bond_count_2_corner2 == 1) or \
               (bond_count_1_corner2 == 1 and bond_count_2_corner2 == 2):
                is_c5_c6 = True
                break
          if is_c5_c6:
              break
        if is_c5_c6:
           has_delta5_bond = True
           break

    if not has_delta5_bond:
        return False, "No double bond between C5 and C6"
    

    return True, "3beta-hydroxy-Delta(5)-steroid identified"