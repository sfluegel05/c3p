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

    # 1. Check for glycosidic bond (C-O-C connecting two anomeric carbons)
    glycosidic_bond_pattern = Chem.MolFromSmarts("[C;H1](-[OX2])-[OX2]-[C;H1](-[OX2])")
    glycosidic_bonds = mol.GetSubstructMatches(glycosidic_bond_pattern)
    if not glycosidic_bonds:
       return False, "No glycosidic bond found"
    if len(glycosidic_bonds) > 1:
        return False, "More than one glycosidic bond found"


    # 2. Check for monosaccharide rings (pyranose and furanose) including the hydroxyl groups
    pyranose_pattern = Chem.MolFromSmarts("[C;H1,H2](-[OX2])-[C;H1,H2](-[OX2])-[C;H1,H2](-[OX2])-[C;H1,H2](-[OX2])-[C;H1,H2](-[OX2])-[O;H0]")
    furanose_pattern = Chem.MolFromSmarts("[C;H1,H2](-[OX2])-[C;H1,H2](-[OX2])-[C;H1,H2](-[OX2])-[C;H1,H2](-[OX2])-[O;H0]")
    hydroxyl_pattern = Chem.MolFromSmarts("[O;H1]-[C]")

    pyranose_matches = mol.GetSubstructMatches(pyranose_pattern)
    furanose_matches = mol.GetSubstructMatches(furanose_pattern)

    if len(pyranose_matches) + len(furanose_matches) != 2:
       return False, "Incorrect number of monosaccharide rings found"
    
    #3. Count total rings
    ring_pattern = Chem.MolFromSmarts("[C;!H0]1-[C;!H0]-[C;!H0]-[C;!H0]-[C;!H0]-1")
    rings = mol.GetSubstructMatches(ring_pattern)
    if len(rings) > 10:
        return False, "Too many rings, likely not a disaccharide"

    # 4. Check if the glycosidic bond connects the two monosaccharide units.
    
    ring_atoms = set()
    for match in pyranose_matches:
        ring_atoms.update(match)
    for match in furanose_matches:
         ring_atoms.update(match)

    atom1_idx = glycosidic_bonds[0][0]
    atom2_idx = glycosidic_bonds[0][2]
    
    
    #Check if both atoms of glycosidic bond are in a ring:
    if not((atom1_idx in ring_atoms) and (atom2_idx in ring_atoms)):
         return False, "Glycosidic bond not connecting two monosaccharide rings"
    
    #Check if they belong to two *different* rings
    
    found_connection=False
    for ring_match in pyranose_matches:
        if (atom1_idx in ring_match and atom2_idx not in ring_match) or (atom2_idx in ring_match and atom1_idx not in ring_match):
          found_connection = True
          break

    if found_connection:
      return True, "Contains two monosaccharide units linked by a glycosidic bond"

    for ring_match in furanose_matches:
        if (atom1_idx in ring_match and atom2_idx not in ring_match) or (atom2_idx in ring_match and atom1_idx not in ring_match):
            found_connection = True
            break
    if not found_connection:
        return False, "Glycosidic bond not connecting two different monosaccharide rings"
    
    return True, "Contains two monosaccharide units linked by a glycosidic bond"