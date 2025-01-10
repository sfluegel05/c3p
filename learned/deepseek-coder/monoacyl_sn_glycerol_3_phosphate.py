"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
"""
Classifies: CHEBI:17408 monoacyl-sn-glycerol 3-phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate based on its SMILES string.
    A monoacyl-sn-glycerol 3-phosphate has a glycerol backbone with a phosphate group at the 3rd position
    and a single acyl group at either the 1st or 2nd position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacyl-sn-glycerol 3-phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone pattern (C-C-C with hydroxyl groups)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for phosphate group attached to the 3rd carbon (more flexible pattern)
    phosphate_pattern = Chem.MolFromSmarts("[CH2X4][OX2][PX4](=[OX1])([OX2H0-1])[OX2H0-1]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found at the 3rd position"

    # Check for a single ester group (acyl group attached via ester bond)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Check for fatty acid chain (long carbon chain attached to ester)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 1:
        return False, "Missing fatty acid chain"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Chain too short to be a fatty acid"

    # Check molecular weight - monoacyl-sn-glycerol 3-phosphate typically >300 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for monoacyl-sn-glycerol 3-phosphate"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 10:
        return False, "Too few carbons for monoacyl-sn-glycerol 3-phosphate"
    if o_count < 5 or o_count > 7:
        return False, "Incorrect number of oxygens for monoacyl-sn-glycerol 3-phosphate"

    return True, "Contains glycerol backbone with a single acyl group and a phosphate group at the 3rd position"