"""
Classifies: CHEBI:17408 monoacylglycerol
"""
"""
Classifies: monoacylglycerol
A glyceride in which any one of the R groups (position not specified) is an acyl group 
while the remaining two R groups can be either H or alkyl groups.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string.
    A monoacylglycerol has a glycerol backbone with one fatty acid chain attached via an ester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for basic glycerol backbone pattern - more flexible version
    glycerol_pattern = Chem.MolFromSmarts("[OX2H1,OX2R]-[CH2X4]-[CHX4]-[CH2X4]-[OX2H1,OX2R]")
    glycerol_pattern2 = Chem.MolFromSmarts("[CH2X4](-[OX2H1,OX2R])-[CHX4](-[OX2H1,OX2R])-[CH2X4](-[OX2H1,OX2R])")
    
    has_glycerol = mol.HasSubstructMatch(glycerol_pattern) or mol.HasSubstructMatch(glycerol_pattern2)
    if not has_glycerol:
        return False, "No glycerol backbone found"

    # Look for exactly one ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[#6]")  # Ensure carbon chain after ester
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Count free hydroxyl groups - should have exactly 2
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) != 2:
        return False, f"Found {len(hydroxyl_matches)} free hydroxyl groups, need exactly 2"

    # Verify the molecule is not primarily a ring system
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 1:  # Allow one ring for some natural products
        return False, "Too many ring systems for monoacylglycerol"

    # Basic element counts and checks
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 5:  # 3 from glycerol + at least 2 from acyl group
        return False, "Too few carbons for monoacylglycerol"
    if o_count != 4:  # Must have exactly 4 oxygens (2 OH + 2 ester)
        return False, "Must have exactly 4 oxygens"

    # Verify proper connectivity
    # The glycerol carbon with the ester must connect to exactly one ester oxygen
    ester_o = mol.GetSubstructMatch(Chem.MolFromSmarts("[OX2][CX3](=[OX1])"))[0]
    for atom in mol.GetAtomWithIdx(ester_o).GetNeighbors():
        if atom.GetAtomicNum() == 6:  # Carbon
            if len([n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]) > 2:
                return False, "Invalid connectivity around glycerol carbon"

    return True, "Contains glycerol backbone with one acyl group attached via ester bond"