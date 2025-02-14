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

    # Define SMARTS pattern for a generic monosaccharide ring (5 or 6 membered with at least one oxygen)
    monosaccharide_ring_smarts = "[CX4;R5,R6]-[CX4;R5,R6]-[CX4;R5,R6]-[CX4;R5,R6]-[CX4;R5,R6]-[OX2;R5,R6]"  # 6-membered
    monosaccharide_ring_smarts_5 = "[CX4;R5]-[CX4;R5]-[CX4;R5]-[CX4;R5]-[OX2;R5]"  # 5-membered
    monosaccharide_pattern = Chem.MolFromSmarts(monosaccharide_ring_smarts)
    monosaccharide_pattern_5 = Chem.MolFromSmarts(monosaccharide_ring_smarts_5)

    if monosaccharide_pattern is None or monosaccharide_pattern_5 is None:
        return None, "Invalid SMARTS pattern, check syntax"

    # Find all matches of the monosaccharide ring pattern
    matches = mol.GetSubstructMatches(monosaccharide_pattern)
    matches_5 = mol.GetSubstructMatches(monosaccharide_pattern_5)
    num_rings = len(matches) + len(matches_5)

    if num_rings == 0:
      return False, "No monosaccharide rings found."

    # Count glycosidic bonds (C-O-C linkages) between rings. The carbons MUST be part of a ring
    glycosidic_bond_pattern = Chem.MolFromSmarts("[CX4;R5,R6]-[OX2]-[CX4;R5,R6]")  #C-O-C linkage where Carbons are part of a ring
    if glycosidic_bond_pattern is None:
        return None, "Invalid SMARTS pattern for glycosidic bonds"
    
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)
    num_glycosidic = len(glycosidic_matches)

    if num_glycosidic == 0 and num_rings > 1:
      return False, "Multiple monosaccharide rings detected, but no glycosidic bonds"

    # Check if there are more than 10 monosaccharide residues, if not, move to other checks
    if num_rings > 10:
      return True, f"Contains {num_rings} monosaccharide residues, linked by glycosidic bonds. Likely a polysaccharide."

    # Check molecular weight - polysaccharides typically > 500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
      return False, "Molecular weight too low for polysaccharide"

    # Heuristic approach:
    if num_rings > 3 and num_glycosidic > 1 :
      return True, f"Multiple monosaccharide rings, linked by glycosidic bonds. Likely an oligosaccharide/polysaccharide."

    # If all checks passed, but still below 10
    if num_rings > 2 and num_glycosidic > 1:
      return True, f"Multiple monosaccharide rings, possibly an oligosaccharide, needs further analysis."

    return False, f"Fails polysaccharide criteria. {num_rings} rings found. Needs further analysis."