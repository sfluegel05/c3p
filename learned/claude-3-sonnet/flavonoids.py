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
    
    # Check carbon count (should be around C15-C16 core + possible substituents)
    if c_count < 15:
        return False, "Too few carbons for flavonoid core structure"
        
    # Must contain oxygen
    if o_count < 2:
        return False, "Insufficient oxygen atoms for flavonoid"

    # Common flavonoid core patterns
    patterns = [
        # Basic flavone core
        "[#6]1~[#6]~[#6](~[#6]~[#6]~[#6]1)~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]~[#6]~[#6]2",
        # Chalcone pattern
        "[#6]~[#6](=O)~[#6]=?[#6]~[#6]1~[#6]~[#6]~[#6]~[#6]~[#6]1",
        # Isoflavone core
        "O=C1C=C(c2ccccc2)Oc2ccccc12",
        # Flavanone core
        "O=C1CC(c2ccccc2)Oc2ccccc12",
        # Aurone core
        "O=C1C=C(c2ccccc2)c2ccccc2O1"
    ]
    
    # Check for presence of any core patterns
    found_core = False
    for pattern in patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_core = True
            break
            
    if not found_core:
        return False, "No flavonoid core structure found"

    # Check for aromatic rings (should have at least 2)
    ring_info = mol.GetRingInfo()
    aromatic_rings = 0
    for ring in ring_info.AtomRings():
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            aromatic_rings += 1
    
    if aromatic_rings < 2:
        return False, "Insufficient aromatic rings for flavonoid structure"

    # Look for typical oxygen-containing functional groups
    o_patterns = [
        # Ketone/carbonyl
        "[#6]=O",
        # Hydroxyl
        "[#6]-[OH]",
        # Ether
        "[#6]-O-[#6]",
        # Ester
        "[#6](=O)-O-[#6]"
    ]
    
    oxygen_groups = 0
    for pattern in o_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            oxygen_groups += 1
            
    if oxygen_groups < 2:
        return False, "Insufficient oxygen-containing functional groups"

    # Additional checks for specific flavonoid characteristics
    # Check for presence of phenyl rings connected by a propane chain or similar
    phenyl_propane = "[#6]1~[#6]~[#6]~[#6]~[#6]~[#6]1-[#6]-[#6]-[#6]"
    if not mol.HasSubstructMatch(Chem.MolFromSmarts(phenyl_propane)):
        # Some flavonoids might not match this exact pattern due to modifications
        # but we'll still consider them if they passed other checks
        pass

    # Calculate number of rings
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 2:
        return False, "Insufficient ring count for flavonoid structure"

    return True, "Contains flavonoid core structure with appropriate substitution patterns"