"""
Classifies: CHEBI:39434 limonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    Limonoids are highly oxygenated triterpenoids with a 4,4,8-trimethyl-17-furanylsteroid skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a limonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for basic molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350 or mol_wt > 850:  # Increased upper limit
        return False, "Molecular weight outside typical limonoid range (350-850)"

    # Count rings - limonoids have multiple fused rings
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    if ring_count < 4:
        return False, "Too few rings for a limonoid (minimum 4)"

    # Look for furan or modified furan-like structures with more variations
    furan_patterns = [
        'c1ccoc1',  # classical furan
        'C1=COC=C1',  # furan alternative representation
        'O1C=CC=C1',  # another furan representation
        'C1=COC(=O)C1',  # furanone
        'O1C=CC(=O)C1',  # isofuranone
        'O1C=C[CH2]C1',  # dihydrofuran
        'O1C(=O)C=CC1',  # furan-2-one
        'c1cocc1'  # another furan representation
    ]
    
    has_furan = False
    for pattern in furan_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_furan = True
            break

    # Core structure patterns - more flexible to catch variations
    core_patterns = [
        # Basic steroid-like core with various possible substitutions
        '[C;R]1~[C;R]~[C;R]~[C;R]2~[C;R]([C;R]1)~[C;R]1~[C;R]~[C;R]~[C;R]~[C;R]~[C;R]12',
        # Alternative pattern with oxygen bridges
        '[C;R]1~[C;R]~[C;R]2~O~[C;R]~[C;R]~[C;R]2~[C;R]1',
        # Pattern for typical A/B ring fusion
        '[C;R]1~[C;R]2~[C;R]~[C;R]~[C;R](~[C;R]1)~[C;R]2',
        # Pattern for B/C ring fusion
        '[C;R]1~[C;R]2~[C;R]~[C;R]~[C;R]~[C;R]2~[C;R]~[C;R]1'
    ]
    
    has_core = False
    for pattern in core_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_core = True
            break

    # Check for oxygenation pattern
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 5:  # Increased minimum oxygen count
        return False, "Insufficient oxygenation for a limonoid"

    # Look for characteristic functional groups with more variations
    functional_groups = [
        ('[C;R]-C(=O)-O-[C,H]', 'ester'),
        ('[C;R]-C(=O)-[C;R]', 'ketone'),
        ('[C;R]-[OH]', 'hydroxyl'),
        ('C(C)(C)([OH])', 'tertiary alcohol'),
        ('O1[C;R][C;R]1', 'epoxide'),
        ('[C;R]-O-C(=O)', 'ester linkage'),
        ('[C;R]-O-[C;R]', 'ether bridge'),
        ('C(=O)-O-[C;R]', 'ester group')
    ]
    
    found_groups = []
    for pattern, group_type in functional_groups:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_groups.append(group_type)
    
    # Count carbons and check for reasonable range for limonoids
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20 or c_count > 40:
        return False, "Carbon count outside typical limonoid range (20-40)"

    # Final decision based on combined criteria
    if has_core and len(found_groups) >= 2 and (has_furan or o_count >= 7):
        return True, f"Matches limonoid characteristics: contains {', '.join(found_groups)}, appropriate ring system" + \
               (" and furan moiety" if has_furan else " and high oxygenation pattern")
    
    if not has_core:
        return False, "Missing characteristic limonoid core structure"
    if len(found_groups) < 2:
        return False, "Insufficient characteristic functional groups"
    if not has_furan and o_count < 7:
        return False, "Missing both furan moiety and sufficient oxygenation"
    
    return False, "Does not match overall limonoid structural requirements"