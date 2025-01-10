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
    if c_count < 15:
        return False, "Too few carbons for flavonoid core structure"
        
    # Must contain oxygen
    if o_count < 2:
        return False, "Insufficient oxygen atoms for flavonoid"

    # Count aromatic rings
    ring_info = mol.GetRingInfo()
    aromatic_rings = 0
    for ring in ring_info.AtomRings():
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            aromatic_rings += 1
    
    if aromatic_rings < 2:
        return False, "Insufficient aromatic rings for flavonoid structure"

    # Core structure patterns for different flavonoid classes
    patterns = [
        # Basic flavone/flavonol core (with variations in saturation)
        "O=C1CC(c2ccccc2)Oc2ccccc12",
        # Chalcone pattern (with variations in saturation)
        "O=C(CC)c1ccccc1",
        # Isoflavone core variation
        "O=C1C=Cc2ccccc2O1",
        # Flavanone core variation
        "O1CCCCC1",
        # General diphenylpropane skeleton
        "c1ccccc1CCc1ccccc1"
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
        "C(=O)",
        # Hydroxyl
        "OH",
        # Ether
        "COC",
        # Pyranone ring (common in flavonoids)
        "O=C1CCOc2ccccc12"
    ]
    
    oxygen_groups = 0
    for pattern in o_patterns:
        pattern_mol = Chem.MolFromSmarts(pattern)
        if pattern_mol is not None and mol.HasSubstructMatch(pattern_mol):
            oxygen_groups += 1
            
    if oxygen_groups < 2:
        return False, "Insufficient oxygen-containing functional groups"

    # Calculate molecular properties
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 200 or mol_weight > 1000:
        return False, "Molecular weight outside typical flavonoid range"

    # Calculate number of rings
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 2:
        return False, "Insufficient ring count for flavonoid structure"

    # Check for conjugated system
    conjugated_pattern = Chem.MolFromSmarts("c1ccccc1")
    if conjugated_pattern is not None and not mol.HasSubstructMatch(conjugated_pattern):
        return False, "Missing conjugated aromatic system"

    # Additional structural features common in flavonoids
    features = [
        # Phenol group
        "c1ccccc1O",
        # Common substitution pattern
        "Oc1ccccc1",
        # Pyrone-type ring
        "O=C1CCOc2ccccc12"
    ]
    
    feature_count = 0
    for pattern in features:
        pattern_mol = Chem.MolFromSmarts(pattern)
        if pattern_mol is not None and mol.HasSubstructMatch(pattern_mol):
            feature_count += 1

    if feature_count == 0:
        return False, "Missing typical flavonoid substitution patterns"

    return True, "Contains flavonoid core structure with appropriate substitution patterns"