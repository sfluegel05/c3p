"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    Glycosaminoglycans are polysaccharides containing a substantial proportion of aminomonosaccharide residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosaminoglycan, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # 1. Check for the presence of multiple 5- or 6-membered rings with oxygen (potential monosaccharide units)
    ring_pattern_5 = Chem.MolFromSmarts("C1OCC[C,O]1") # 5-membered ring
    ring_pattern_6 = Chem.MolFromSmarts("C1COCCC[C,O]1") # 6-membered ring
    ring_matches_5 = mol.GetSubstructMatches(ring_pattern_5)
    ring_matches_6 = mol.GetSubstructMatches(ring_pattern_6)
    total_rings = len(ring_matches_5) + len(ring_matches_6)
    if total_rings < 3:
       return False, f"Too few monosaccharide rings: {total_rings}"

    # 2. Check for glycosidic linkages (C-O-C) connecting rings (at least 2)
    glycosidic_link_pattern = Chem.MolFromSmarts("[CX4]-[OX2]-[CX4]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_link_pattern)
    if len(glycosidic_matches) < 2:
        return False, f"Too few glycosidic linkages: {len(glycosidic_matches)}"

    # 3. Check for presence of amino groups (-NH2) or substituted amino groups (-NHR).
    amino_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(N-[CX4]=O)]") # NH2 or NHR, excluding amide nitrogens
    amino_matches = mol.GetSubstructMatches(amino_pattern)
    if len(amino_matches) < 1:
      #Check for N-acetylated amino groups (common modification)
      n_acetyl_pattern = Chem.MolFromSmarts("[NX3;!H0]-[CX3](=[OX1])-[CX4]")
      n_acetyl_matches = mol.GetSubstructMatches(n_acetyl_pattern)
      if len(n_acetyl_matches) < 1:
        return False, f"No amino groups or N-acetylated amino groups detected."
      else:
         amino_count = len(n_acetyl_matches)
    else:
      amino_count = len(amino_matches)

    #4.  Check for enough carbons.  GAG's are large molecules with many carbon atoms.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 12:
      return False, "Too few carbons for a glycosaminoglycan"


    # 5. Check for a minimum number of oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 5:
        return False, "Too few oxygens for glycosaminoglycan"


    # 6. Check for minimum number of hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 3:
        return False, "Too few hydroxyl groups for a glycosaminoglycan"


    if total_rings >= 3 and len(glycosidic_matches) >= 2 and amino_count >= 1:
       return True, "Contains multiple monosaccharide units, glycosidic bonds, and amino groups."
    else:
       return None, None