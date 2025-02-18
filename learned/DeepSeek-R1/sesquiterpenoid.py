"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    Sesquiterpenoids are C15-based terpenoids that may have structural modifications.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Basic carbon count check with wider range for modifications
    num_c = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_c < 12 or num_c > 30:  # Expanded upper limit for functional groups
        return False, f"Carbon count {num_c} outside 12-30 range"
    
    # Check for characteristic methyl groups (minimum 2)
    methyl_pattern = Chem.MolFromSmarts("[CH3;!$(C=O)]")  # Exclude carbonyl-attached methyls
    methyl_matches = len(mol.GetSubstructMatches(methyl_pattern))
    if methyl_matches < 2:
        return False, f"Only {methyl_matches} methyl groups"
    
    # Terpenoid structural features check
    ring_info = mol.GetRingInfo()
    rings = ring_info.NumRings()
    double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    
    # Require either ring systems or multiple double bonds typical of terpenoids
    if rings < 1 and double_bonds < 2:
        return False, "Lacks terpenoid features (rings/double bonds)"
    
    # Check for terpene-like branching patterns using more flexible smarts
    # Looks for 3 consecutive carbons with at least one branch (common in terpenes)
    branch_pattern = Chem.MolFromSmarts("[CH2][CH]([!H])[CH2]")
    if not mol.HasSubstructMatch(branch_pattern):
        return False, "No characteristic branching pattern"
    
    return True, f"{num_c} carbons, {methyl_matches} methyl groups, terpenoid features"