"""
Classifies: CHEBI:48953 cyclohexenones
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdSubstructLibrary


def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule is a cyclohexenone based on its SMILES string.
    A cyclohexenone is a six-membered alicyclic ketone having one double bond in the ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclohexenone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns
    six_membered_ring = Chem.MolFromSmarts("C1CCCCC1")
    ketone_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    double_bond_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")

    # Check for 6 membered ring
    if not mol.HasSubstructMatch(six_membered_ring):
      return False, "Molecule does not contain a six-membered ring."

    # Check for ketone and double bond in the ring
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    
    if not ketone_matches or not double_bond_matches:
        return False, "Molecule does not contain ketone or double bond groups"
    
    ring_atoms = mol.GetSubstructMatch(six_membered_ring)
    ketone_in_ring = False
    for match in ketone_matches:
        for atom_index in match:
            if atom_index in ring_atoms:
                ketone_in_ring = True
                break
        if ketone_in_ring:
            break
    if not ketone_in_ring:
        return False, "Ketone is not in ring"

    double_bond_in_ring = False
    for match in double_bond_matches:
      for atom_index in match:
        if atom_index in ring_atoms:
            double_bond_in_ring = True
            break
      if double_bond_in_ring:
          break
    if not double_bond_in_ring:
        return False, "Double bond is not in ring"

    
    # Count double bonds in the ring.
    double_bond_count = 0
    for match in double_bond_matches:
        
        atom_in_ring = False
        for atom_index in match:
            if atom_index in ring_atoms:
                atom_in_ring = True
                break
        
        if atom_in_ring:
           double_bond_count += 1

    if double_bond_count != 1:
        return False, f"Must have exactly one double bond in ring, got {double_bond_count}"
    
    return True, "Molecule is a cyclohexenone"