"""
Classifies: CHEBI:4194 D-hexose
"""
"""
Classifies: CHEBI:16234 D-hexose
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    D-hexose is a hexose (6-carbon monosaccharide) with D-configuration at C5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a D-hexose, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic composition check
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count != 6:
        return False, f"Must have exactly 6 carbons, found {c_count}"
    if o_count != 6:
        return False, f"Must have exactly 6 oxygens, found {o_count}"

    # Check for non C/H/O atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [1, 6, 8]:
            return False, "Contains elements other than C, H, O"

    # Check for carboxyl groups (would indicate uronic acid)
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=O)[OX2H1,OX1-]")):
        return False, "Contains carboxyl group"

    # Check for basic monosaccharide structure
    # Pattern for open chain form
    aldehyde_pattern = Chem.MolFromSmarts("[CH1](=O)[CH1][CH1][CH1][CH1][CH2]")
    # Pattern for pyranose form (6-membered ring)
    pyranose_pattern = Chem.MolFromSmarts("[CH2]1[CH1][CH1][CH1][CH1][OH1]1")
    # Pattern for furanose form (5-membered ring)
    furanose_pattern = Chem.MolFromSmarts("[CH2]1[CH1][CH1][CH1][OH1]1")
    
    is_aldehyde = mol.HasSubstructMatch(aldehyde_pattern)
    is_pyranose = mol.HasSubstructMatch(pyranose_pattern)
    is_furanose = mol.HasSubstructMatch(furanose_pattern)
    
    if not (is_aldehyde or is_pyranose or is_furanose):
        return False, "Does not match basic monosaccharide structure"

    # Check for D-configuration at C5 position
    # For pyranose forms
    d_pyranose_pattern = Chem.MolFromSmarts("[CH2][C@H]1O[CH1][CH1][CH1][CH1]1")
    # For furanose forms
    d_furanose_pattern = Chem.MolFromSmarts("[CH2][C@H]1O[CH1][CH1][CH1]1")
    # For open chain forms
    d_aldehyde_pattern = Chem.MolFromSmarts("[CH2][C@H](O)[CH1](O)[CH1](O)[CH1](O)C=O")
    
    # Count chiral centers
    chiral_centers = Chem.FindMolChiralCenters(mol)
    if len(chiral_centers) < 4:
        return False, "Insufficient chiral centers for a D-hexose"

    # Check for D-configuration patterns
    if mol.HasSubstructMatch(d_pyranose_pattern) or \
       mol.HasSubstructMatch(d_furanose_pattern) or \
       mol.HasSubstructMatch(d_aldehyde_pattern):
        
        # Additional check for hydroxyl groups
        hydroxyls = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[OH1]")))
        if hydroxyls < 5:
            return False, "Insufficient hydroxyl groups"
            
        if is_pyranose:
            return True, "D-hexose in pyranose form"
        elif is_furanose:
            return True, "D-hexose in furanose form"
        else:
            return True, "D-hexose in open chain form"
            
    return False, "Does not have D-configuration at C5"