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

    # Basic molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350 or mol_wt > 850:
        return False, "Molecular weight outside typical limonoid range (350-850)"

    # More specific limonoid core patterns
    core_patterns = [
        # Basic 4-ring steroid-like core with methyl positions
        '[C;R]1[C;R][C;R]2[C;R](C)[C;R][C;R][C@]2(C)[C;R]2[C;R][C;R][C;R]3[C;R][C;R][C;R][C;R]3[C;R]12',
        # Alternative pattern with oxygen bridge
        '[C;R]1[C;R]2O[C;R][C;R]([C;R]2)[C;R]2[C;R][C;R][C;R]3[C;R][C;R][C;R][C;R]3[C;R]12',
        # Pattern for modified A/B/C/D ring system
        '[C;R]1[C;R]2[C;R][C;R]3[C;R][C;R][C;R]4[C;R][C;R][C;R][C;R]4[C;R]3[C;R]2[C;R][C;R]1'
    ]
    
    has_core = False
    for pattern in core_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_core = True
            break

    if not has_core:
        return False, "Missing characteristic limonoid core structure"

    # Check for furan ring specifically attached to the core
    furan_patterns = [
        # Furan connected to specific position
        '[C;R]1[C;R][C;R]2[C;R][C;R][C;R][C;R]2[C;R]2[C;R][C;R](c3ccoc3)[C;R][C;R]12',
        # Alternative furan connection
        '[C;R]([C;R]1[C;R][C;R]2)([C;R]2[C;R])c3ccoc3',
        # Modified furan pattern
        '[C;R]1[C;R][C;R]2[C;R][C;R][C;R][C;R]2[C;R]2[C;R][C;R](C3=COC=C3)[C;R][C;R]12'
    ]
    
    has_furan = False
    for pattern in furan_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_furan = True
            break

    # Check for characteristic oxygen-containing groups
    o_patterns = {
        'lactone': '[C;R]1[C;R]C(=O)O[C;R]1',
        'epoxy': '[C;R]1O[C;R]1',
        'ketone': '[C;R]C(=O)[C;R]',
        'acetoxy': 'CC(=O)O[C;R]',
        'hydroxy': '[C;R][OH]',
        'ether bridge': '[C;R]O[C;R]'
    }
    
    found_groups = []
    for group, pattern in o_patterns.items():
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_groups.append(group)

    # Count oxygens but don't use as strict criterion
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Check for characteristic methyl groups
    methyl_pattern = Chem.MolFromSmarts('[C;R]([C;R])(C)(C)')
    has_gem_dimethyl = mol.HasSubstructMatch(methyl_pattern)

    # Decision making with more weight on structural features
    if has_core and has_furan and len(found_groups) >= 2:
        return True, f"Matches limonoid characteristics: contains furan ring and {', '.join(found_groups)}"
    
    if has_core and len(found_groups) >= 3 and o_count >= 6:
        return True, f"Matches limonoid characteristics: highly oxygenated with {', '.join(found_groups)}"
        
    if not has_furan and o_count < 5:
        return False, "Insufficient oxygenation and missing furan ring"
        
    if len(found_groups) < 2:
        return False, "Insufficient characteristic functional groups"

    return False, "Does not match overall limonoid structural requirements"