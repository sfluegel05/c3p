"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
"""
Classifies: Triterpenoid saponins
Definition: A terpene glycoside in which the terpene moiety is a triterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a triterpenoid saponin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbons to check if it's a triterpenoid (should have ~30 carbons in core)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 25:  # Allow some flexibility
        return False, f"Too few carbons ({carbon_count}) for a triterpenoid core"
    
    # Look for glycosidic linkage patterns - multiple SMARTS to catch different cases
    glycosidic_patterns = [
        Chem.MolFromSmarts("[CR]-[OR2]-[CR]"),  # Basic glycosidic bond
        Chem.MolFromSmarts("[CR]-O-[C;R1][C;R1][C;R1]"),  # Connection to sugar ring
        Chem.MolFromSmarts("[CR]-O-[C;R1]1O[C;R1][C;R1][C;R1][C;R1]1")  # More specific sugar pattern
    ]
    has_glycosidic = False
    for pattern in glycosidic_patterns:
        if mol.HasSubstructMatch(pattern):
            has_glycosidic = True
            break
    if not has_glycosidic:
        return False, "No glycosidic linkage found"
    
    # Look for sugar moiety patterns (multiple patterns to catch variations)
    sugar_patterns = [
        Chem.MolFromSmarts("[CR1]1O[CR1][CR1][CR1][CR1]1"), # Basic pyranose
        Chem.MolFromSmarts("[CR1]1O[CR1][CR1][CR1][CR1][CR1]1"), # 6-membered sugar
        Chem.MolFromSmarts("[CR1]1O[CR1]([CR1][CR1][CR1]1)O") # Sugar with OH group
    ]
    has_sugar = False
    for pattern in sugar_patterns:
        if mol.HasSubstructMatch(pattern):
            has_sugar = True
            break
    if not has_sugar:
        return False, "No sugar moiety found"
    
    # Look for triterpenoid core patterns (multiple connected rings)
    triterpenoid_patterns = [
        Chem.MolFromSmarts("C1~C~C~C~C~1~C1~C~C~C~C~1"), # Basic two rings
        Chem.MolFromSmarts("C12CCC3C(C1)CCC4C3(C)CCC2C4"), # More specific steroid-like core
        Chem.MolFromSmarts("C1CC2CCC3(C)C(C2)C(CC4C3(C)CCC4)C1") # Alternative core pattern
    ]
    has_triterpenoid = False
    for pattern in triterpenoid_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            has_triterpenoid = True
            break
    if not has_triterpenoid:
        return False, "Missing characteristic triterpenoid ring pattern"
    
    # Check ring count (triterpenoids typically have 4-6 rings)
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    if ring_count < 4:
        return False, f"Too few rings ({ring_count}) for a triterpenoid"
    
    # Count oxygens (should have multiple due to glycosylation and hydroxyl groups)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 5:  # Reduced threshold
        return False, f"Too few oxygens ({oxygen_count}) for a triterpenoid saponin"
    
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:  # Lowered threshold
        return False, f"Molecular weight ({mol_wt:.1f}) too low for triterpenoid saponin"
    
    # Count hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_matches < 2:
        return False, f"Too few hydroxyl groups ({hydroxyl_matches}) for a saponin"
        
    return True, "Contains triterpenoid core with glycosidic linkages and characteristic structural features"