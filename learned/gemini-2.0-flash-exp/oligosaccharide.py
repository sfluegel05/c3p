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
      
    # Define a more general SMARTS pattern for monosaccharide rings
    # Cyclic structure of carbons and one oxygen, with at least 3 hydroxy groups
    sugar_ring_pattern = Chem.MolFromSmarts("[C;R]1([OH])[C;R]([OH])[C;R]([OH])[C;R][O;R]1")

    # Count potential sugar rings
    sugar_ring_matches = mol.GetSubstructMatches(sugar_ring_pattern)

    if len(sugar_ring_matches) < 2:
      return False, "Molecule does not contain at least two sugar rings."

    # Verify that rings are connected by glycosidic linkages and are carbohydrates
    
    verified_rings = 0
    for match in glycosidic_matches:
        start_atom = mol.GetAtomWithIdx(match[0])
        mid_atom = mol.GetAtomWithIdx(match[1])
        end_atom = mol.GetAtomWithIdx(match[2])

        if not (mid_atom.GetAtomicNum()==8 and start_atom.IsInRing() and end_atom.IsInRing()):
          continue
        
        #verify if this is a real carbohydrate link by checking for -OH groups
        is_valid_start = False
        is_valid_end = False
        
        for ring_match in sugar_ring_matches:
          for ring_atom_idx in ring_match:
            if mol.GetAtomWithIdx(ring_atom_idx).IsInRing() and ring_atom_idx == start_atom.GetIdx():
              is_valid_start = True
            if mol.GetAtomWithIdx(ring_atom_idx).IsInRing() and ring_atom_idx == end_atom.GetIdx():
              is_valid_end = True
        
        if is_valid_start and is_valid_end:
          verified_rings += 2
        
    if verified_rings < 4:
        return False, "Not all rings are properly connected by glycosidic linkages."

    # Check molecular weight - oligosaccharides typically >300 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for oligosaccharide"

    return True, "Molecule contains multiple sugar rings linked by glycosidic bonds."