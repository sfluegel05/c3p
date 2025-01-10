"""
Classifies: CHEBI:16749 1-phosphatidyl-1D-myo-inositol
"""
"""
Classifies: CHEBI:28874 1-phosphatidyl-1D-myo-inositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_phosphatidyl_1D_myo_inositol(smiles: str):
    """
    Determines if a molecule is a 1-phosphatidyl-1D-myo-inositol based on its SMILES string.
    A 1-phosphatidyl-1D-myo-inositol consists of a glycerol backbone with two fatty acid chains,
    a phosphate group, and a 1D-myo-inositol moiety attached to the phosphate group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-phosphatidyl-1D-myo-inositol, False otherwise
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
        
    # Look for 2 ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"

    # Look for phosphate group (-O-P(=O)(-O)-)
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2])[OX2]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Look for 1D-myo-inositol moiety (6-membered ring with 6 hydroxyl groups)
    inositol_pattern = Chem.MolFromSmarts("[C@H]1([C@@H]([C@H]([C@@H]([C@H]([C@@H]1O)O)O)O)O)")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No 1D-myo-inositol moiety found"

    # Check that the phosphate group is attached to the inositol moiety
    phosphate_inositol_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2])[OX2][C@H]1([C@@H]([C@H]([C@@H]([C@H]([C@@H]1O)O)O)O)O)")
    if not mol.HasSubstructMatch(phosphate_inositol_pattern):
        return False, "Phosphate group not attached to 1D-myo-inositol moiety"

    return True, "Contains glycerol backbone with two fatty acid chains, a phosphate group, and a 1D-myo-inositol moiety attached to the phosphate group"