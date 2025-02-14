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
        return False, "Invalid SMILES string"

    # 1. Define a pattern for a monosaccharide ring (5 or 6 membered) with a hemiacetal/ketal.
    monosaccharide_pattern = Chem.MolFromSmarts("[C1][O][C]([C])([C])[C][C]1") # 5 or 6 membered ring.  This pattern captures the hemiacetal carbon and its neighbours (simplification)
    monosaccharide_matches = mol.GetSubstructMatches(monosaccharide_pattern)
    
    #If no monosaccharide pattern is detected, it cannot be a GAG
    if len(monosaccharide_matches) == 0:
        return False, "No monosaccharide ring detected"

    # 2. Check for amino groups (NH2, NHR, or N-acetylated) *on or adjacent to* the monosaccharide rings
    amino_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(N-[CX4]=O)]") #NH2 or NHR, excluding amide nitrogens
    n_acetyl_pattern = Chem.MolFromSmarts("[NX3;!H0]-[CX3](=[OX1])-[CX4]")# N-acetyl group

    amino_matches = []
    for match in monosaccharide_matches:
       for atom_index in match: # check atoms in each monosaccharide
            monosaccharide_atom = mol.GetAtomWithIdx(atom_index)
            for neighbor in monosaccharide_atom.GetNeighbors():
                if neighbor.Match(amino_pattern) or neighbor.Match(n_acetyl_pattern):
                    amino_matches.append(neighbor.GetIdx())

    if len(amino_matches) == 0:
        return False, f"No amino groups or N-acetylated amino groups detected on or adjacent to monosaccharide ring(s)."

    # 3. Check for at least one glycosidic linkage (C-O-C) connected to at least one of the monosaccharide units.
    glycosidic_link_pattern = Chem.MolFromSmarts("[CX4]-[OX2]-[CX4]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_link_pattern)
    
    has_glycosidic = False
    for match in glycosidic_matches:
        for atom_index in match:
          glycosidic_atom = mol.GetAtomWithIdx(atom_index)
          for monosaccharide_match in monosaccharide_matches:
              for mono_atom_index in monosaccharide_match:
                  mono_atom = mol.GetAtomWithIdx(mono_atom_index)
                  if glycosidic_atom.IsInRing() and (glycosidic_atom.GetIdx() in [x.GetIdx() for x in mono_atom.GetNeighbors()]):
                    has_glycosidic = True
                    break
              if has_glycosidic:
                  break
        if has_glycosidic:
            break;
    if not has_glycosidic:
       return False, "No glycosidic linkage detected, connected to a monosaccharide unit"

    return True, "Contains monosaccharide units, glycosidic bonds, and amino groups on or adjacent to monosaccharides."