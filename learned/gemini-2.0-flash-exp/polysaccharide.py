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

    # Check molecular weight - polysaccharides typically > 500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
      return False, f"Molecular weight too low for polysaccharide"
    
    # More specific SMARTS pattern for a monosaccharide ring
    #   - Six membered ring with 1 O and 5 C or 5 membered ring with 1 O and 4C
    #   - At least 2 hydroxyl groups
    #   - Check for hemiacetal or acetal oxygen
    monosaccharide_ring_smarts = "[CX4H2,CX4H1,CX4;R1]([OX2;R1])[CX4H2,CX4H1;R1]([CX4H2,CX4H1;R1])([CX4H2,CX4H1;R1])[CX4H2,CX4H1;R1][OX2H;R1]" # 6-membered ring with one oxygen and at least one OH
    monosaccharide_ring_smarts2 = "[CX3H1,CX3;R1]([OX2;R1])[CX3H1;R1]([CX3H1;R1])([CX3H1;R1])[OX2H;R1]" # 5-membered ring with one oxygen and at least one OH
    monosaccharide_pattern1 = Chem.MolFromSmarts(monosaccharide_ring_smarts)
    monosaccharide_pattern2 = Chem.MolFromSmarts(monosaccharide_ring_smarts2)
    if monosaccharide_pattern1 is None or monosaccharide_pattern2 is None:
        return None, "Invalid SMARTS pattern, check syntax"

    # Find all matches of the monosaccharide ring patterns
    matches1 = mol.GetSubstructMatches(monosaccharide_pattern1)
    matches2 = mol.GetSubstructMatches(monosaccharide_pattern2)
    num_rings = len(matches1) + len(matches2)


    # Count glycosidic bonds (C-O-C linkages) between monosaccharide rings
    glycosidic_bond_pattern = Chem.MolFromSmarts("[CX4;!H0;R1]([OX2;R1])[CX4;!H0;R2]")
    if glycosidic_bond_pattern is None:
        return None, "Invalid SMARTS pattern for glycosidic bonds"
    
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)
    num_glycosidic = len(glycosidic_matches)


    # A polysaccharide should have at least 10 monosaccharide units
    # we estimate the number of monosaccharide units as the number of glycosidic bonds +1
    num_monosaccharide_units = num_glycosidic + 1
    if num_monosaccharide_units >= 10 or num_rings >= 10:
       return True, f"Contains {num_rings} monosaccharide rings and at least {num_monosaccharide_units} glycosidic bonds. Likely a polysaccharide."
    

    # Heuristic approach for potential polysaccharides
    if num_rings > 3 and num_glycosidic > 1:
      return True, f"Multiple monosaccharide rings, linked by glycosidic bonds, suggesting an oligosaccharide/polysaccharide. {num_rings} rings, {num_glycosidic} bonds"

    # If all checks passed, but still below 10
    if num_rings > 2 and num_glycosidic > 1:
      return True, f"Multiple monosaccharide rings, possibly an oligosaccharide. Needs further analysis. {num_rings} rings, {num_glycosidic} bonds"
      
    if num_rings == 0:
      return False, f"No monosaccharide rings found. {num_rings} rings, {num_glycosidic} bonds"

    if num_glycosidic == 0 and num_rings > 1 :
        return False, f"Multiple monosaccharide rings detected, but no glycosidic bonds. {num_rings} rings, {num_glycosidic} bonds"


    return False, f"Fails polysaccharide criteria. {num_rings} rings, {num_glycosidic} bonds. Needs further analysis."