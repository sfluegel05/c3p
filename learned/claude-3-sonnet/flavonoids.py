"""
Classifies: CHEBI:72544 flavonoids
"""
"""
Classifies: flavonoids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_flavonoid, reason)
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic element counts
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Check carbon count (C15/C16 core + possible substituents)
    if c_count < 12:  # Relaxed from 15 to catch simpler flavonoids
        return False, "Too few carbons for flavonoid core structure"
        
    # Must contain oxygen
    if o_count < 1:  # Relaxed from 2
        return False, "Insufficient oxygen atoms for flavonoid"

    # Core structure patterns for different flavonoid classes
    patterns = [
        # Basic flavone/flavonol core (more flexible)
        "[#6]1~[#6]~[#6](~[#6]2~[#6]~[#6]~[#6]~[#6]~[#6]2)~[#8]~[#6]2~[#6]~[#6]~[#6]~[#6]~[#6]12",
        # Isoflavone core
        "[#6]1~[#6]~[#6]2~[#8]~[#6]~[#6](=O)~[#6]2=C~[#6]1",
        # Chalcone pattern
        "[#6]1~[#6]~[#6]~[#6]~[#6]~[#6]1~[#6](=O)~[#6]~[#6]~[#6]1~[#6]~[#6]~[#6]~[#6]~[#6]1",
        # Flavanone core
        "[#6]1~[#6]~[#6](~[#6]2~[#6]~[#6]~[#6]~[#6]~[#6]2)~[#8]~[#6]2~[#6](=O)~[#6]~[#6]~[#6]12",
        # Anthocyanidin core
        "[#6]1~[#6]~[#6]2~[#8]~[#6](~[#6]3~[#6]~[#6]~[#6]~[#6]~[#6]3)~[#6](~[#8])~[#6]~[#6]2~[#6]~[#6]1",
        # Aurone core
        "[#6]1~[#6]~[#6]2~[#8]~[#6](=C~[#6]3~[#6]~[#6]~[#6]~[#6]~[#6]3)~[#6](=O)~[#6]2~[#6]~[#6]1",
        # More general diphenylpropane skeleton
        "[#6]1~[#6]~[#6]~[#6]~[#6]~[#6]1~[#6]~[#6]~[#6]1~[#6]~[#6]~[#6]~[#6]~[#6]1"
    ]
    
    # Check for presence of any core patterns
    found_core = False
    for pattern in patterns:
        pattern_mol = Chem.MolFromSmarts(pattern)
        if pattern_mol is not None and mol.HasSubstructMatch(pattern_mol):
            found_core = True
            break
            
    if not found_core:
        return False, "No flavonoid core structure found"

    # Look for typical oxygen-containing functional groups
    o_patterns = [
        # Ketone/carbonyl
        "[#6][#6](=[#8])[#6]",
        # Hydroxyl
        "[#8H1]",
        # Ether
        "[#6][#8][#6]",
        # Glycoside pattern
        "[#6][#8][#6]1[#8][#6]([#6][#8])[#6]([#8])[#6]([#8])[#6]1",
        # Pyranone ring
        "[#6]1[#6][#6](=[#8])[#8][#6]2[#6][#6][#6][#6][#6]12"
    ]
    
    oxygen_groups = 0
    for pattern in o_patterns:
        pattern_mol = Chem.MolFromSmarts(pattern)
        if pattern_mol is not None:
            matches = mol.GetSubstructMatches(pattern_mol)
            oxygen_groups += len(matches)
            
    if oxygen_groups == 0:  # Relaxed from 2
        return False, "No characteristic oxygen-containing functional groups"

    # Additional structural features common in flavonoids
    features = [
        # Phenol group
        "c1ccccc1[OH]",
        # Common substitution patterns
        "[#8]c1ccccc1",
        # Conjugated system
        "c1ccccc1",
        # Common glycosylation pattern
        "[#6][#8][#6]1[#8][#6][#6][#6][#6][#6]1[#8]"
    ]
    
    feature_count = 0
    for pattern in features:
        pattern_mol = Chem.MolFromSmarts(pattern)
        if pattern_mol is not None and mol.HasSubstructMatch(pattern_mol):
            feature_count += 1

    if feature_count == 0:
        return False, "Missing typical flavonoid structural features"

    # Calculate number of rings
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 2:
        return False, "Insufficient ring count for flavonoid structure"

    return True, "Contains flavonoid core structure with appropriate substitution patterns"