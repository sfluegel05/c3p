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

    # 1D-myo-inositol ring pattern with specific stereochemistry
    # The pattern matches the characteristic 1D-myo-inositol with:
    # - axial OH at position 2
    # - equatorial OHs at positions 3,4,5,6
    # - position 1 connected to phosphate
    myo_inositol_pattern = Chem.MolFromSmarts(
        "[O,P]-[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)1"
    )
    
    if not mol.HasSubstructMatch(myo_inositol_pattern):
        return False, "No 1D-myo-inositol ring with correct stereochemistry found"

    # Phosphate group connected to position 1 of inositol
    phosphate_pattern = Chem.MolFromSmarts(
        "[O,OH]-[P](=O)([O,OH])-O-[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)1"
    )
    
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group at position 1 of inositol"

    # Glycerophospholipid backbone with two ester-linked fatty acids
    glycerol_pattern = Chem.MolFromSmarts(
        "[#6]C(=O)OC[CH](COC(=O)[#6])COP(=O)([O,OH])-O"
    )
    
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerophospholipid backbone found"

    # Check for two ester groups
    ester_pattern = Chem.MolFromSmarts("[#6]C(=O)O[#6]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"

    # Check for fatty acid chains (at least 4 carbons in length)
    fatty_acid_pattern = Chem.MolFromSmarts("[#6]~[#6]~[#6]~[#6]")
    fatty_acid_matches = len(mol.GetSubstructMatches(fatty_acid_pattern))
    if fatty_acid_matches < 2:
        return False, "Missing fatty acid chains"

    # Verify overall composition
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    if p_count != 1:
        return False, "Must have exactly one phosphorus atom"
    if o_count < 12:  # 6 from inositol, 2 from esters, 4 from phosphate
        return False, "Insufficient oxygen atoms"
    if c_count < 20:  # minimum for shortest possible chains
        return False, "Carbon chains too short"

    return True, "Contains 1D-myo-inositol ring with phosphate at position 1 and two fatty acid chains"