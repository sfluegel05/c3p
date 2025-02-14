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

    # 1. Check for glycosidic bonds (C-O-C connecting two anomeric carbons)
    glycosidic_bond_pattern = Chem.MolFromSmarts("[C;H1,H2](-[OX2])-[OX2]-[C;H1,H2](-[OX2])")
    glycosidic_bonds = mol.GetSubstructMatches(glycosidic_bond_pattern)

    if not glycosidic_bonds:
        return False, "No glycosidic bond found"

    # 2. Check for two monosaccharide rings
    pyranose_pattern = Chem.MolFromSmarts("[C;H1,H2](-[OX2])-[C;H1,H2](-[OX2])-[C;H1,H2](-[OX2])-[C;H1,H2](-[OX2])-[C;H1,H2](-[OX2])-O")
    furanose_pattern = Chem.MolFromSmarts("[C;H1,H2](-[OX2])-[C;H1,H2](-[OX2])-[C;H1,H2](-[OX2])-[C;H1,H2](-[OX2])-O")
    
    pyranose_matches = mol.GetSubstructMatches(pyranose_pattern)
    furanose_matches = mol.GetSubstructMatches(furanose_pattern)

    if len(pyranose_matches) + len(furanose_matches) < 2:
        return False, "Less than two monosaccharide rings found"

    # 3. Check if the glycosidic bond connect the two monosaccharide units. 
    # Here we're using the indexes from matches, and checking if a glycosidic bond is between two different sugar rings
    
    ring_atoms = set()
    for match in pyranose_matches:
        ring_atoms.update(match)
    for match in furanose_matches:
        ring_atoms.update(match)
    
    
    valid_glycosidic_bonds = 0

    for match in glycosidic_bonds:
        atom1_idx = match[0]
        atom2_idx = match[2]
        if (atom1_idx in ring_atoms) and (atom2_idx in ring_atoms):
                
            #Check if they belong to two *different* rings
            for ring_match in pyranose_matches:
               if (atom1_idx in ring_match and atom2_idx not in ring_match) or (atom2_idx in ring_match and atom1_idx not in ring_match):
                   valid_glycosidic_bonds += 1
                   break
            
            if valid_glycosidic_bonds > 0:
              continue
            for ring_match in furanose_matches:
               if (atom1_idx in ring_match and atom2_idx not in ring_match) or (atom2_idx in ring_match and atom1_idx not in ring_match):
                   valid_glycosidic_bonds += 1
                   break

    if valid_glycosidic_bonds == 0:
        return False, "No glycosidic bond connecting two monosaccharide rings found."

    return True, "Contains two monosaccharide units linked by a glycosidic bond"