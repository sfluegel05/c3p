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
    if carbon_count < 20:  # Lowered threshold
        return False, f"Too few carbons ({carbon_count}) for a triterpenoid core"
    
    # Define common triterpenoid core patterns with more flexible matching
    triterpenoid_patterns = [
        # Basic pentacyclic pattern (more flexible)
        "*~1~*~*~2~*~3~*~*~4~*~5~*~*~4~*~*~5~*~*~3~*~*~2~*~*~1",
        # Tetracyclic pattern
        "*~1~*~*~2~*~3~*~*~4~*~*~3~*~*~4~*~*~2~*~*~1",
        # Oleanane/Ursane core (more flexible)
        "*~1~*~*~2~*~3~*(*~*~4~*~5~*~*~4~*)~*~*~3~*~*~2~*~*~1",
        # Modified core with lactone
        "*~1~*~*~2~*~3~*~*~4~*~5~*~*~4~*~*~5~*~*~3~*~*~2~O~*~1",
        # Dammarane type
        "*~1~*~*~2~*~3~*~*~4~*~*~3~*~*~4~*~*~2~*~*~1",
        # More general fused ring system
        "*~1~*~*~2~*~3~*~*~4~*~*~3~*~*~4~*~*~2~*~*~1"
    ]
    
    has_triterpenoid = False
    for pattern in triterpenoid_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt is not None and mol.HasSubstructMatch(patt):
            has_triterpenoid = True
            break
    
    if not has_triterpenoid:
        return False, "No characteristic triterpenoid ring system found"
    
    # Look for sugar patterns (more flexible)
    sugar_patterns = [
        # General pyranose pattern
        "O1[C;R1][C;R1][C;R1][C;R1][C;R1]1",
        # Furanose pattern
        "O1[C;R1][C;R1][C;R1][C;R1]1",
        # More flexible sugar pattern
        "O1[C;R1][C;R1]([O,C])[C;R1]([O,C])[C;R1]1",
        # Pattern for modified sugars
        "O1[C;R1][C;R1]([O,C])[C;R1]([O,C])[C;R1]([O,C])O1"
    ]
    
    has_sugar = False
    for pattern in sugar_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt is not None and mol.HasSubstructMatch(patt):
            has_sugar = True
            break
            
    if not has_sugar:
        return False, "No sugar moiety found"
    
    # Look for glycosidic linkage (more flexible patterns)
    glycosidic_patterns = [
        # General glycosidic linkage
        "[C;R0,R1]-[O;R0]-[C;R1]",
        # Specific O-glycosidic bond
        "[C;R0,R1]-[O;R0]-[C;R1]1[O;R1][C;R1][C;R1][C;R1][C;R1]1",
        # Acetal linkage
        "[C;R1]1[O;R1][C;R1][C;R1][C;R1][C;R1]1[O;R0][C;R0,R1]"
    ]
    
    has_glycosidic = False
    for pattern in glycosidic_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt is not None and mol.HasSubstructMatch(patt):
            has_glycosidic = True
            break
            
    if not has_glycosidic:
        return False, "No glycosidic linkage found"
    
    # Check for characteristic functional groups
    functional_groups = {
        "hydroxyl": "[OH]",
        "carboxyl": "[C;!R](=O)[OH]",
        "ketone": "[C;!R,R](=O)[C;!R,R]",
        "methyl": "[CH3]",
        "ether": "[C;!R,R]-[O;R0]-[C;!R,R]"
    }
    
    # Count functional groups
    group_counts = {}
    for name, pattern in functional_groups.items():
        patt = Chem.MolFromSmarts(pattern)
        if patt is not None:
            group_counts[name] = len(mol.GetSubstructMatches(patt))
    
    # Check ring count
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    if ring_count < 4:
        return False, f"Too few rings ({ring_count}) for a triterpenoid"
    
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:  # Lowered threshold
        return False, f"Molecular weight ({mol_wt:.1f}) too low for triterpenoid saponin"
    
    # Verify minimum requirements for functional groups
    if group_counts.get("hydroxyl", 0) < 2:
        return False, "Too few hydroxyl groups for a saponin"
    
    if group_counts.get("methyl", 0) < 3:
        return False, "Too few methyl groups for a triterpenoid"
        
    return True, "Contains triterpenoid core with glycosidic linkages and characteristic structural features"