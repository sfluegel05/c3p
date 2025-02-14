"""
Classifies: CHEBI:50699 oligosaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Get number of rings
    ring_info = mol.GetRingInfo()
    n_rings = ring_info.NumRings()

    if n_rings < 2:
        return False, "Molecule has less than 2 rings, not an oligosaccharide"
    
    # Check for glycosidic linkages (C-O-C between rings)
    glycosidic_linkage_pattern = Chem.MolFromSmarts("[CX4][OX2][CX4]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_linkage_pattern)
    
    ring_atoms_count = 0
    for atom in mol.GetAtoms():
       if atom.IsInRing():
         ring_atoms_count+=1
         if atom.GetAtomicNum() != 6 and atom.GetAtomicNum() != 8:
            return False, "Rings must contain only carbon and oxygen atoms."


    if len(glycosidic_matches) < 1:
        return False, "Molecule does not contain any glycosidic linkages."

    # Check that the molecule consists mainly of 5- and 6-membered rings with multiple -OH groups.
    sugar_ring_pattern1 = Chem.MolFromSmarts("[C1]([C][C][C][C][O1])[O]") #5 member ring
    sugar_ring_pattern2 = Chem.MolFromSmarts("[C1]([C][C][C][C][C][O1])[O]") #6 member ring
    
    sugar_matches1 = mol.GetSubstructMatches(sugar_ring_pattern1)
    sugar_matches2 = mol.GetSubstructMatches(sugar_ring_pattern2)
    if len(sugar_matches1) + len(sugar_matches2) == 0:
         return False, "Molecule does not contain sugar rings (5 or 6 membered cycles with an oxygen in the ring)."

    #Check number of carbons and oxygen
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if c_count < 6: #minimum 2 monosaccharides are required, hence 2*3 carbons
       return False, "Oligosaccharide should have at least 6 carbon atoms"

    if o_count < 4: #minimum 2 monosaccharides are required, hence 2*2 oxygens
       return False, "Oligosaccharide should have at least 4 oxygen atoms"
    
    # Check the linkages between cycles
    for match in glycosidic_matches:
      start_atom = mol.GetAtomWithIdx(match[0])
      mid_atom = mol.GetAtomWithIdx(match[1])
      end_atom = mol.GetAtomWithIdx(match[2])
      if not (mid_atom.GetAtomicNum()==8 and start_atom.IsInRing() and end_atom.IsInRing()):
         return False, "Linkages between the cycles must be ether links."
         

    return True, "Molecule contains multiple sugar rings linked by glycosidic bonds."