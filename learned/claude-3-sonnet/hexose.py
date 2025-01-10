"""
Classifies: CHEBI:18133 hexose
"""
"""
Classifies: CHEBI:24895 hexose
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    A hexose is a six-carbon monosaccharide with either an aldehyde group (aldohexose)
    or a ketone group (ketohexose).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count != 6:
        return False, f"Must have exactly 6 carbons, found {carbon_count}"

    # Count oxygens (hexoses typically have 6 oxygens in cyclic form, 5 in linear form)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 5 or oxygen_count > 6:
        return False, f"Expected 5-6 oxygens, found {oxygen_count}"

    # Check for hydroxyl groups (at least 4 needed)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    if mol.HasSubstructMatch(hydroxyl_pattern):
        hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
        if hydroxyl_count < 4:
            return False, f"Too few hydroxyl groups ({hydroxyl_count}), minimum 4 required"
    else:
        return False, "No hydroxyl groups found"

    # Check for common sugar ring patterns (furanose or pyranose)
    furanose_pattern = Chem.MolFromSmarts("[C]1[C][C][C](O1)")
    pyranose_pattern = Chem.MolFromSmarts("[C]1[C][C][C][C](O1)")
    
    # Check for aldehyde or ketone group
    aldehyde_pattern = Chem.MolFromSmarts("[CH1](=O)")
    ketone_pattern = Chem.MolFromSmarts("[CX3](=O)")
    
    has_ring = mol.HasSubstructMatch(furanose_pattern) or mol.HasSubstructMatch(pyranose_pattern)
    has_carbonyl = mol.HasSubstructMatch(aldehyde_pattern) or mol.HasSubstructMatch(ketone_pattern)
    
    if not (has_ring or has_carbonyl):
        return False, "No sugar ring or carbonyl group found"

    # Check it's not a polysaccharide by looking for glycosidic bonds
    glycosidic_pattern = Chem.MolFromSmarts("[OX2]([CH1][OX2])[CH1]")
    if mol.HasSubstructMatch(glycosidic_pattern):
        glycosidic_count = len(mol.GetSubstructMatches(glycosidic_pattern))
        if glycosidic_count > 1:
            return False, "Appears to be a polysaccharide"

    # Additional check for ester groups (would indicate a modified sugar)
    ester_pattern = Chem.MolFromSmarts("[#6][CX3](=O)[OX2][#6]")
    if mol.HasSubstructMatch(ester_pattern):
        return False, "Contains ester group(s) - modified sugar"

    # Calculate ring count
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 1:
        return False, "Contains multiple rings"

    return True, "Matches hexose pattern with correct number of carbons and hydroxyls"