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
    
    # Look for glycosidic linkage pattern (C-O-C where one C is part of a ring)
    glycosidic_pattern = Chem.MolFromSmarts("[CR1,CR2,CR3]-[OR2]-[CR1,CR2,CR3;R1]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    if not glycosidic_matches:
        return False, "No glycosidic linkage found"
    
    # Look for sugar moiety patterns (pyranose rings)
    sugar_pattern = Chem.MolFromSmarts("[CR1][OR2][CR1][CR1][CR1][CR1][OR2]")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No sugar moiety found"
    
    # Check for multiple rings (triterpenoids typically have 4-6 rings)
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    if ring_count < 4:
        return False, f"Too few rings ({ring_count}) for a triterpenoid"
    
    # Check for multiple hydroxyl groups (characteristic of saponins)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, "Too few hydroxyl groups for a saponin"
    
    # Look for characteristic triterpenoid patterns (multiple connected rings)
    triterpenoid_pattern = Chem.MolFromSmarts("C1~C~C~C~C~1~C1~C~C~C~C~1")  # Simplified ring pattern
    if not mol.HasSubstructMatch(triterpenoid_pattern):
        return False, "Missing characteristic triterpenoid ring pattern"
    
    # Calculate molecular weight (should be relatively high due to glycosylation)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:  # Typical triterpenoid saponins are >600 Da
        return False, f"Molecular weight ({mol_wt:.1f}) too low for triterpenoid saponin"
    
    # Count oxygens (should have multiple due to glycosylation and hydroxyl groups)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 6:
        return False, f"Too few oxygens ({oxygen_count}) for a triterpenoid saponin"
    
    # Calculate number of rotatable bonds (should be relatively high due to sugar moieties)
    rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if rotatable_bonds < 5:
        return False, "Too few rotatable bonds for glycosylated structure"
        
    return True, "Contains triterpenoid core with glycosidic linkages and characteristic structural features"