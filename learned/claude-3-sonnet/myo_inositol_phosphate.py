"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
"""
Classifies: myo-inositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a myo-inositol phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for cyclohexane core
    cyclohexane_pattern = Chem.MolFromSmarts("[C@H]1[C@H][C@H][C@H][C@H][C@H]1")
    if not mol.HasSubstructMatch(cyclohexane_pattern):
        return False, "No cyclohexane core found"
    
    # Check for phosphate groups (-OP(=O)(O)O or -OP(=O)(O)[O-])
    phosphate_pattern = Chem.MolFromSmarts("[OX2][P](=[OX1])([OX2,OX1-])[OX2,OX1-]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) == 0:
        return False, "No phosphate groups found"
    
    # Count carbons (should be exactly 6 for myo-inositol)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 6:
        return False, f"Found {c_count} carbons, need exactly 6 for myo-inositol"
    
    # Count oxygens (should be at least 6 + 4 per phosphate)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    min_o_count = 6 + (4 * len(phosphate_matches))  # 6 for hydroxyl/phosphate attachment + 4 per phosphate
    if o_count < min_o_count:
        return False, f"Insufficient oxygen atoms for structure"
    
    # Check for phosphorus
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count == 0:
        return False, "No phosphorus atoms found"
    if p_count != len(phosphate_matches):
        return False, "Inconsistent phosphate group count"
        
    # Check for ring size (should be 6)
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    if not any(size == 6 for size in ring_sizes):
        return False, "No 6-membered ring found"
        
    # Check for substitution pattern
    # All carbons should have one oxygen (either OH or OP)
    carbon_oxygen_pattern = Chem.MolFromSmarts("[#6]~[OX2]")
    carbon_oxygen_matches = mol.GetSubstructMatches(carbon_oxygen_pattern)
    if len(carbon_oxygen_matches) < 6:
        return False, "Not all carbons have oxygen substituents"
    
    return True, f"Contains myo-inositol core with {len(phosphate_matches)} phosphate groups"