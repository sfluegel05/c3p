"""
Classifies: CHEBI:50699 oligosaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors


def is_oligosaccharide(smiles: str):
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string.
    An oligosaccharide is a compound formed by a few monosaccharide units joined by glycosidic linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oligosaccharide, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycosidic linkages (C-O-C between rings)
    glycosidic_linkage_pattern = Chem.MolFromSmarts("[C;R][O;R][C;R]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_linkage_pattern)
    if not glycosidic_matches:
      return False, "Molecule does not contain any glycosidic linkages (C-O-C between rings)."

    # Define SMARTS patterns for common monosaccharide rings
    sugar_ring_pattern_5 = Chem.MolFromSmarts("[C;R1,R2]1[C;R1,R2]([OH])[C;R1,R2]([OH])[C;R1,R2]([OH])[O;R1,R2]1") # 5 membered ring
    sugar_ring_pattern_6 = Chem.MolFromSmarts("[C;R1,R2]1[C;R1,R2]([OH])[C;R1,R2]([OH])[C;R1,R2]([OH])[C;R1,R2]([OH])[O;R1,R2]1") # 6 membered ring

    #Get all rings
    ring_info = mol.GetRingInfo()
    n_rings = ring_info.NumRings()
    if n_rings < 2:
        return False, "Molecule has less than 2 rings, not an oligosaccharide"

    # Count sugar rings
    sugar_ring_matches_5 = mol.GetSubstructMatches(sugar_ring_pattern_5)
    sugar_ring_matches_6 = mol.GetSubstructMatches(sugar_ring_pattern_6)
    n_sugar_rings = len(sugar_ring_matches_5) + len(sugar_ring_matches_6)
    
    if n_sugar_rings == 0:
        return False, "Molecule does not contain any sugar rings."

    # Check if most rings are sugar rings (allow for a few non-sugar rings, in case of modifications)
    if n_sugar_rings < n_rings * 0.7 :
       return False, "Molecule does not have a majority of sugar rings."

    # Check all ring atoms are either C or O
    for atom in mol.GetAtoms():
      if atom.IsInRing() and atom.GetAtomicNum() != 6 and atom.GetAtomicNum() != 8:
         return False, "Rings must contain only carbon and oxygen atoms."

    #Check the link between rings
    for match in glycosidic_matches:
        start_atom = mol.GetAtomWithIdx(match[0])
        mid_atom = mol.GetAtomWithIdx(match[1])
        end_atom = mol.GetAtomWithIdx(match[2])
        if not (mid_atom.GetAtomicNum()==8 and start_atom.IsInRing() and end_atom.IsInRing()):
            return False, "Linkages between the cycles must be ether links between two rings."

    # Check molecular weight - oligosaccharides typically >300 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for oligosaccharide"


    return True, "Molecule contains multiple sugar rings linked by glycosidic bonds."