"""
Classifies: CHEBI:36233 disaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide consists of two monosaccharides linked by a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a disaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for glycosidic bond (C-O-C)
    glycosidic_bond_pattern = Chem.MolFromSmarts("[CX4]-[OX2]-[CX4]")
    glycosidic_bonds = mol.GetSubstructMatches(glycosidic_bond_pattern)
    
    if not glycosidic_bonds:
       return False, "No glycosidic bond found"
    if len(glycosidic_bonds) > 1:
        return False, "More than one glycosidic bond found"
    
    # Get atoms involved in the glycosidic bond
    atom1_idx = glycosidic_bonds[0][0]
    atom2_idx = glycosidic_bonds[0][2]

    # 2. Find monosaccharide rings (5 or 6 membered)
    ring_pattern = Chem.MolFromSmarts("[C;!H0]1-[C;!H0]-[C;!H0]-[C;!H0]-[C;!H0]-1") # 5-membered ring
    ring_pattern2 = Chem.MolFromSmarts("[C;!H0]1-[C;!H0]-[C;!H0]-[C;!H0]-[C;!H0]-[C;!H0]-1") # 6-membered ring

    rings = mol.GetSubstructMatches(ring_pattern)
    rings2 = mol.GetSubstructMatches(ring_pattern2)
    
    if len(rings) + len(rings2) != 2:
       return False, "Incorrect number of monosaccharide rings found"
    
    # Collect ring atoms
    ring_atoms = []
    for match in rings:
        ring_atoms.append(set(match))
    for match in rings2:
        ring_atoms.append(set(match))
        
    #Check if both atoms of glycosidic bond are in a ring:
    found_connection = False
    for ring_set in ring_atoms:
         if (atom1_idx in ring_set) and (atom2_idx in ring_set):
            return False, "Glycosidic bond connecting atoms in the same ring"
         if (atom1_idx in ring_set) or (atom2_idx in ring_set):
           found_connection = True
    
    if not found_connection:
          return False, "Glycosidic bond not connecting any monosaccharide ring"

    #Check if they belong to two *different* rings
    
    found_connection = False
    
    for i in range(len(ring_atoms)):
        for j in range(i+1, len(ring_atoms)):
           if (atom1_idx in ring_atoms[i] and atom2_idx in ring_atoms[j]) or (atom2_idx in ring_atoms[i] and atom1_idx in ring_atoms[j]):
                found_connection=True
                break
        if found_connection:
            break
    if not found_connection:
         return False, "Glycosidic bond not connecting two different monosaccharide rings"
   
    # Check for anomeric carbons (carbon attached to two oxygens)
    anomeric_pattern = Chem.MolFromSmarts("C([OX2])-[OX2]")
    anomeric_matches = mol.GetSubstructMatches(anomeric_pattern)

    # Check that the carbons of the glycosidic bonds are anomeric
    anomeric_carbons = set()
    for match in anomeric_matches:
         anomeric_carbons.update(match)

    if not (atom1_idx in anomeric_carbons and atom2_idx in anomeric_carbons):
        return False, "Glycosidic bond not connecting two anomeric carbons"

    return True, "Contains two monosaccharide units linked by a glycosidic bond"