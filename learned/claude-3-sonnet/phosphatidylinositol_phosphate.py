"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
"""
Classifies: phosphatidylinositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate based on its SMILES string.
    PIPs contain a myo-inositol ring with at least one phosphate group attached directly
    to the ring. Most common forms are PI(3)P, PI(4)P, PI(5)P, PI(3,4)P2, PI(3,5)P2,
    PI(4,5)P2, and PI(3,4,5)P3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a PIP, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for inositol ring with at least one OH/OP group
    inositol_pattern = Chem.MolFromSmarts("[CH1]1[CH1][CH1][CH1][CH1][CH1]1")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"

    # Look for phosphate groups attached to inositol ring
    # Pattern matches inositol carbon with attached phosphate
    inositol_phosphate_pattern = Chem.MolFromSmarts("[C]1[C][C][C][C][C]1[OX2]P(=O)([OX2])[OX2]")
    phosphate_on_ring = mol.GetSubstructMatches(inositol_phosphate_pattern)
    
    if not phosphate_on_ring:
        return False, "No phosphate groups attached to inositol ring"

    # Count total phosphorus atoms
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    # Count oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if p_count < 1:
        return False, "No phosphorus atoms found"
    
    if o_count < 8:
        return False, "Too few oxygen atoms for a PIP"

    # Look for either:
    # 1. Standard glycerol-based PIP structure
    # 2. Modified PIP structure (like the hexakis ester)
    
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = len(mol.GetSubstructMatches(ester_pattern))
    
    # Check for standard glycerol-based structure
    if mol.HasSubstructMatch(glycerol_pattern):
        # Verify fatty acid chains for glycerol-based PIPs
        fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
        fatty_acid_matches = len(mol.GetSubstructMatches(fatty_acid_pattern))
        
        if fatty_acid_matches < 2 and ester_matches < 2:
            return False, "Missing fatty acid chains in glycerol-based PIP"
    
    # For non-glycerol PIPs, check for sufficient modifications
    elif ester_matches < 1:
        return False, "No ester groups found in modified PIP"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:  # Lowered threshold to accommodate modified PIPs
        return False, "Molecular weight too low for PIP"

    # Classify based on number of phosphates
    if p_count >= 4:
        return True, "Contains inositol ring with trisphosphate modification (PIP3)"
    elif p_count >= 3:
        return True, "Contains inositol ring with bisphosphate modification (PIP2)"
    else:
        return True, "Contains inositol ring with monophosphate modification (PIP)"