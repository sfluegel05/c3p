"""
Classifies: CHEBI:17408 monoacylglycerol
"""
"""
Classifies: monoacylglycerol
A glyceride with one acyl group attached to a glycerol backbone
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

    # Look for glycerol backbone pattern (C-C-C with oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[OX2H1,OX2R]-[CH2X4]-[CHX4]-[CH2X4]-[OX2H1,OX2R]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
        
    # Look for exactly one ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Count hydroxyl groups - should have exactly 2 free OH groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) != 2:
        return False, f"Found {len(hydroxyl_matches)} free hydroxyl groups, need exactly 2"

    # Basic element counts
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Minimum requirements for a monoacylglycerol
    if c_count < 5:  # 3 from glycerol + at least 2 from acyl group
        return False, "Too few carbons for monoacylglycerol"
    if o_count < 4:  # 2 OH groups + 2 from ester
        return False, "Too few oxygens for monoacylglycerol"

    # Check for reasonable molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 120:  # Approximate minimum weight for smallest monoacylglycerol
        return False, "Molecular weight too low for monoacylglycerol"

    # Verify carbon chain attached to ester
    acyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[CX4]")
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No carbon chain attached to ester group"

    return True, "Contains glycerol backbone with one acyl group attached via ester bond"