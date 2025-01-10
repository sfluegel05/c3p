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
    if carbon_count < 25:
        return False, f"Too few carbons ({carbon_count}) for a triterpenoid core"
    
    # Define common triterpenoid core patterns
    triterpenoid_patterns = [
        # Oleanane type core (pentacyclic)
        "[C]1[C][C]2[C]3[C]([C]4[C]([C]5[C]([C][C]4)[C]([C][C]5)[C])[C][C]3)[C][C]2[C][C]1",
        # Ursane type core
        "[C]1[C][C]2[C]3[C]([C]4[C]([C]5[C]([C][C]4)[C]([C][C]5)[C])[C][C]3)[C][C]2[C][C]1",
        # Dammarane type core (tetracyclic)
        "[C]1[C][C]2[C]3[C]([C]4[C]([C][C]3)[C]([C][C]4)[C])[C][C]2[C][C]1",
        # Cycloartane type core
        "[C]1[C][C]2[C]3[C]([C]4[C]([C]5[C]([C][C]4)[C]([C][C]5))[C][C]3)[C][C]2[C][C]1",
        # More general fused ring pattern for triterpenoids
        "[C]1[C][C]2[C]3[C]([C]4[C]([C][C]3)[C]([C][C]4))[C][C]2[C][C]1"
    ]
    
    has_triterpenoid = False
    for pattern in triterpenoid_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt is not None and mol.HasSubstructMatch(patt):
            has_triterpenoid = True
            break
    
    if not has_triterpenoid:
        return False, "No characteristic triterpenoid ring system found"
    
    # Look for sugar patterns with specific stereochemistry
    sugar_patterns = [
        # Pyranose sugar with hydroxyl groups
        "O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)O1",
        # More general pyranose pattern
        "O1[C][C]([O,C])[C]([O,C])[C]([O,C])[C]1[O,C]",
        # Furanose sugar pattern
        "O1[C][C]([O,C])[C]([O,C])[C]1[O,C]"
    ]
    
    has_sugar = False
    for pattern in sugar_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt is not None and mol.HasSubstructMatch(patt):
            has_sugar = True
            break
            
    if not has_sugar:
        return False, "No sugar moiety found"
    
    # Look for glycosidic linkage
    glycosidic_patterns = [
        # O-glycosidic bond
        "[C;R0]-[O;R0]-[C;R1]1[O;R1][C;R1][C;R1][C;R1][C;R1]1",
        # More general glycosidic bond pattern
        "[C;R]-[O;R0]-[C;R1]"
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
        "ketone": "[C;!R](=O)[C;!R]",
        "methyl": "[CH3]"
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
    if mol_wt < 600:  # Increased threshold for triterpenoid saponins
        return False, f"Molecular weight ({mol_wt:.1f}) too low for triterpenoid saponin"
    
    # Verify minimum requirements for functional groups
    if group_counts.get("hydroxyl", 0) < 3:
        return False, "Too few hydroxyl groups for a saponin"
    
    return True, "Contains triterpenoid core with glycosidic linkages and characteristic structural features"