"""
Classifies: CHEBI:16749 1-phosphatidyl-1D-myo-inositol
"""
"""
Classifies: CHEBI:16813 1-phosphatidyl-1D-myo-inositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_phosphatidyl_1D_myo_inositol(smiles: str):
    """
    Determines if a molecule is a 1-phosphatidyl-1D-myo-inositol based on its SMILES string.
    
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

    # Multiple possible SMARTS patterns for myo-inositol ring with correct stereochemistry
    myo_inositol_patterns = [
        # Pattern 1: Direct representation
        "[C@@H]1([OH])[C@H]([OH])[C@@H]([OH])[C@H]([OH])[C@@H]([OH])[C@H]1[OH]",
        # Pattern 2: Alternative representation
        "[C@H]1([OH])[C@@H]([OH])[C@H]([OH])[C@@H]([OH])[C@H]([OH])[C@@H]1[OH]",
        # Pattern 3: More general pattern for cyclohexane with 6 OH groups
        "[CH]1([OH])[CH]([OH])[CH]([OH])[CH]([OH])[CH]([OH])[CH]1[OH]"
    ]
    
    found_inositol = False
    for pattern in myo_inositol_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_inositol = True
            break
            
    if not found_inositol:
        return False, "No myo-inositol ring found or incorrect stereochemistry"

    # Check for phosphate group connected to inositol (-P(=O)(O)-O-)
    phosphate_pattern = Chem.MolFromSmarts("[O][P](=[O])([O,OH])[O]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Check for glycerol backbone with correct stereochemistry
    glycerol_patterns = [
        "[CH2X4][C@H][CH2X4]",
        "[CH2X4][C@@H][CH2X4]"
    ]
    found_glycerol = False
    for pattern in glycerol_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_glycerol = True
            break
            
    if not found_glycerol:
        return False, "No glycerol backbone found or incorrect stereochemistry"

    # Check for two ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])-[#6]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # Check for fatty acid chains (long carbon chains)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, "Missing fatty acid chains"

    # Check connectivity: phosphate must be at position 1 of inositol
    phosphoinositol_pattern = Chem.MolFromSmarts("[O][P](=[O])([O,OH])O[C@H]1[C@H]([OH])[C@@H]([OH])[C@H]([OH])[C@@H]([OH])[C@H]1[OH]")
    if not mol.HasSubstructMatch(phosphoinositol_pattern):
        return False, "Phosphate not correctly connected to position 1 of inositol"

    # Verify overall composition
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    if p_count != 1:
        return False, "Must have exactly one phosphorus atom"
    if o_count < 12:  # 6 from inositol, 2 from esters, 4 from phosphate
        return False, "Insufficient oxygen atoms"
    if c_count < 25:  # minimum for shortest possible chains
        return False, "Carbon chains too short"

    # Check molecular weight - these are typically large molecules
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low"

    return True, "Contains 1D-myo-inositol ring connected via phosphate at position 1 to glycerol with two fatty acid chains"