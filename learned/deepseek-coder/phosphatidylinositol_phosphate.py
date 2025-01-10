"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
"""
Classifies: CHEBI:28874 phosphatidylinositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate based on its SMILES string.
    A phosphatidylinositol phosphate has a glycerol backbone with two fatty acid chains,
    a phosphate group, and a myo-inositol ring with at least one additional phosphate group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylinositol phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 2 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
        
    # Look for 2 ester groups (-O-C(=O)-) for fatty acids
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"

    # Look for phosphate group attached to glycerol
    glycerol_phosphate_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]OP(O)(=O)")
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "No phosphate group attached to glycerol"

    # Look for inositol ring (6 carbons with multiple hydroxyls, allowing for substitutions)
    inositol_pattern = Chem.MolFromSmarts("C1C(O)C(O)C(O)C(O)C(O)C1")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"

    # Look for at least one additional phosphate group on inositol
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(=O)")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 2:
        return False, "Need at least two phosphate groups (one on glycerol and one on inositol)"

    # Check molecular weight - PIPs typically >600 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, "Molecular weight too low for phosphatidylinositol phosphate"

    # Count carbons, oxygens, and phosphates
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    if c_count < 20:
        return False, "Too few carbons for phosphatidylinositol phosphate"
    if o_count < 10:
        return False, "Too few oxygens for phosphatidylinositol phosphate"
    if p_count < 2:
        return False, "Too few phosphates for phosphatidylinositol phosphate"

    return True, "Contains glycerol backbone with two fatty acids, phosphate group, and inositol with additional phosphate"