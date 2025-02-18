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

    # More flexible glycerol backbone patterns that account for stereochemistry
    glycerol_patterns = [
        # Basic glycerol backbone
        "[OX2,OH1]-[CH2X4]-[CHX4]-[CH2X4]-[OX2,OH1]",
        # Pattern matching stereochemistry
        "[OX2,OH1]-[CH2X4]-[C@H,C@@H,CH]-[CH2X4]-[OX2,OH1]",
        # Alternative connection pattern
        "[CH2X4](-[OX2,OH1])-[CHX4](-[OX2,OH1])-[CH2X4](-[OX2,OH1])"
    ]
    
    has_glycerol = False
    for pattern in glycerol_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_glycerol = True
            break
            
    if not has_glycerol:
        return False, "No glycerol backbone found"

    # Look for exactly one ester group (-O-C(=O)-)
    # More flexible ester pattern that accounts for various substitutions
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[#6]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Count substitutable positions (OH or OR groups)
    oh_pattern = Chem.MolFromSmarts("[OX2H1,OX2R]-[CH2X4,CHX4]")
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    if len(oh_matches) < 2:
        return False, "Insufficient hydroxyl or alkoxy groups"

    # Basic element counts and checks
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 5:  # Minimum carbons needed
        return False, "Too few carbons for monoacylglycerol"
    if o_count < 4:  # Minimum oxygens needed
        return False, "Too few oxygens for monoacylglycerol"

    # Check for reasonable molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 120 or mol_wt > 1000:  # Reasonable range for monoacylglycerols
        return False, "Molecular weight outside reasonable range for monoacylglycerol"

    # Verify the acyl chain is attached to the glycerol backbone
    acyl_pattern = Chem.MolFromSmarts("[OX2]-[CH2X4,CHX4]-[CHX4]-[CH2X4]-[OX2][CX3](=[OX1])")
    if not mol.HasSubstructMatch(acyl_pattern):
        # Try alternative pattern
        acyl_pattern2 = Chem.MolFromSmarts("[CX3](=[OX1])[OX2]-[CH2X4,CHX4]-[CHX4](-[OX2,OH1])-[CH2X4][OX2,OH1]")
        if not mol.HasSubstructMatch(acyl_pattern2):
            return False, "Acyl chain not properly connected to glycerol backbone"

    return True, "Contains glycerol backbone with one acyl group attached via ester bond"