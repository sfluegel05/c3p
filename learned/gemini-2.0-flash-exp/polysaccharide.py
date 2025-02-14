"""
Classifies: CHEBI:18154 polysaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polysaccharide(smiles: str):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.
    A polysaccharide is a polymer of monosaccharides linked by glycosidic bonds,
    with more than 10 monosaccharide residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polysaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for a generic monosaccharide ring (5 or 6 membered with at least 2 hydroxyls or oxygens)
    monosaccharide_ring_smarts = "[CX4;R5,R6](-[OX2H1])-[CX4;R5,R6](-[OX2H1])-[CX4;R5,R6]-[CX4;R5,R6]-[CX4;R5,R6]-[OX2;R5,R6]"  # Modified to include ring info
    monosaccharide_pattern = Chem.MolFromSmarts(monosaccharide_ring_smarts)
    if monosaccharide_pattern is None:
        return None, "Invalid SMARTS pattern, check syntax"

    # Find all matches of the monosaccharide ring pattern
    matches = mol.GetSubstructMatches(monosaccharide_pattern)
    num_rings = len(matches)

    if num_rings == 0:
      return False, "No monosaccharide rings found."

    # Count glycosidic bonds (C-O-C linkages) between rings
    glycosidic_bond_pattern = Chem.MolFromSmarts("[CX4;R5,R6]-[OX2]-[CX4;R5,R6]")  #C-O-C linkage between rings
    if glycosidic_bond_pattern is None:
        return None, "Invalid SMARTS pattern for glycosidic bonds"
    
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)
    num_glycosidic = len(glycosidic_matches)

    if num_glycosidic == 0 and num_rings > 1:
      return False, "Multiple monosaccharide rings detected, but no glycosidic bonds"

    # Check the overall carbon and oxygen count
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count < 10 or o_count < 5:
        return False, f"Too few carbons or oxygens to be polysaccharide"

    # Check if there are more than 10 monosaccharide residues
    if num_rings > 10:
      return True, f"Contains {num_rings} monosaccharide residues, linked by glycosidic bonds. Likely a polysaccharide."
    
    # Heuristic approach: count rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if num_rings > 5 and n_rotatable < 5:
      return False, "Too few rotatable bonds for a polysaccharide, despite multiple rings"

    if num_rings > 5 and num_glycosidic < 3:
      return False, "Too few glycosidic bonds compared to rings, not a polysaccharide"
    
    if num_rings > 3 and n_rotatable < 4:
         return False, f"Too few rotatable bonds, needs more than 4."
    
    if num_rings < 4:
      return False, f"Too few monosaccharide rings to be a polysaccharide."

    # If all checks passed, but still below 10
    if num_rings > 2 and num_glycosidic > 1 and n_rotatable > 4 :
      return True, f"Multiple monosaccharide rings, possibly an oligosaccharide, needs further analysis."

    return False, f"Fails polysaccharide criteria. {num_rings} rings found. Needs further analysis."