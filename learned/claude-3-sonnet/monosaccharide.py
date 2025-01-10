"""
Classifies: CHEBI:35381 monosaccharide
"""
"""
Classifies: CHEBI:35381 monosaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a monosaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbons - must have at least 3
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 3:
        return False, "Less than 3 carbon atoms"
        
    # Count oxygens - must have multiple oxygen atoms
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 2:
        return False, "Too few oxygen atoms"
    
    # Look for hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_matches < 2:
        return False, "Insufficient hydroxyl groups"

    # Check for aldehyde group (for aldoses)
    aldehyde_pattern = Chem.MolFromSmarts("[CH,CH2]=O")
    has_aldehyde = mol.HasSubstructMatch(aldehyde_pattern)
    
    # Check for ketone group (for ketoses)
    ketone_pattern = Chem.MolFromSmarts("[#6]-C(=O)-[#6]")
    has_ketone = mol.HasSubstructMatch(ketone_pattern)
    
    # Check for hemiacetal/hemiketal group (for cyclic forms)
    hemiacetal_pattern = Chem.MolFromSmarts("[OH]C1[O,C][C,O][C,O][C,O][C,O]1")
    hemiketal_pattern = Chem.MolFromSmarts("[OH]C1[O,C][C,O][C,O][C,O]1")
    has_cyclic_form = mol.HasSubstructMatch(hemiacetal_pattern) or mol.HasSubstructMatch(hemiketal_pattern)
    
    # Must have either aldehyde/ketone (open form) or hemiacetal/hemiketal (cyclic form)
    if not (has_aldehyde or has_ketone or has_cyclic_form):
        return False, "No aldehyde, ketone, or cyclic sugar form found"
    
    # Check molecular weight - should typically be under 300 Da for single sugar unit
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 300 and not (carbon_count <= 7):  # Allow slightly larger molecules if they have typical sugar carbon count
        return False, "Molecular weight too high for single sugar unit"
    
    # Check for disaccharide linkage pattern
    disaccharide_pattern = Chem.MolFromSmarts("[OX2]([#6])[#6]-[OX2]-[#6]")
    if mol.HasSubstructMatch(disaccharide_pattern):
        # Further check if this is just part of a cyclic structure
        if len(mol.GetSubstructMatches(disaccharide_pattern)) > 1:
            return False, "Contains glycosidic bonds (appears to be oligosaccharide)"
    
    # Calculate ring count
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 1:
        return False, "Contains multiple rings - likely not a monosaccharide"
        
    # Additional check for carbon-hydroxyl ratio (most carbons should have an OH)
    if hydroxyl_matches < (carbon_count - 2) and not has_cyclic_form:
        return False, "Insufficient hydroxyl groups for carbon count"

    return True, "Contains required monosaccharide features (multiple hydroxyls and aldehyde/ketone/cyclic form)"