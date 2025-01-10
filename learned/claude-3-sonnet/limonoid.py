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
    if mol_wt < 350 or mol_wt > 800:
        return False, "Molecular weight outside typical limonoid range (350-800)"

    # Count rings - limonoids have multiple fused rings
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    if ring_count < 4 or ring_count > 8:
        return False, "Number of rings outside typical limonoid range (4-8)"

    # Look for furan or modified furan-like structures
    furan_patterns = [
        'c1ccoc1',  # classical furan
        'C1=COC=C1',  # furan alternative representation
        'O1C=CC=C1',  # another furan representation
        'C1=COC(=O)C1',  # furanone (oxidized furan)
        'O1C=CC(=O)C1'   # isofuranone
    ]
    
    has_furan = False
    for pattern in furan_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_furan = True
            break
    
    if not has_furan:
        return False, "No furan or furan-derived ring found"

    # Check for characteristic limonoid core structure
    # Look for connected 6-membered rings with specific substitution patterns
    core_patterns = [
        # Basic steroid-like tetracyclic core with methyl groups
        '[C;R]1[C;R][C;R][C;R]2[C;R]([C;R]1)(C)[C;R]1[C;R][C;R][C;R]3[C;R](C)[C;R][C;R][C;R]3[C;R]12C',
        # Alternative core pattern with oxygen bridges
        '[C;R]1[C;R][C;R]2O[C;R]3[C;R](O)[C;R][C;R][C;R]3[C;R]2[C;R]1',
        # Pattern for ring D with furan attachment point
        '[C;R]1[C;R][C;R]2[C;R][C;R][C;R](c3ccoc3)[C;R]2[C;R]1'
    ]
    
    has_core = False
    for pattern in core_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_core = True
            break
    
    if not has_core:
        return False, "Missing characteristic limonoid core structure"

    # Check for oxygenation pattern
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 4:
        return False, "Insufficient oxygenation for a limonoid"

    # Look for characteristic functional groups
    functional_groups = [
        ('[C;R]-C(=O)-O-[C;R]', 'ester bridge'),
        ('[C;R]-C(=O)-[C;R]', 'ketone'),
        ('[C;R]-[OH]', 'hydroxyl'),
        ('C(C)(C)([OH])', 'tertiary alcohol'),
        ('O1[C;R][C;R]1', 'epoxide')
    ]
    
    found_groups = []
    for pattern, group_type in functional_groups:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_groups.append(group_type)
    
    if len(found_groups) < 2:
        return False, "Insufficient characteristic functional groups"

    # Count methyl groups in specific positions
    methyl_pattern = Chem.MolFromSmarts('[C;R]-[CH3]')
    methyl_count = len(mol.GetSubstructMatches(methyl_pattern))
    if methyl_count < 3:
        return False, "Insufficient methyl groups for characteristic limonoid substitution"

    return True, f"Matches limonoid characteristics: contains {' and '.join(found_groups)}, appropriate ring system, and furan moiety"