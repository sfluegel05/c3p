from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

def is_abietane_diterpenoid(smiles: str):
    """
    Determines if a molecule is an abietane diterpenoid based on its structure.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an abietane diterpenoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
        
    # Check for C20 diterpenoid core
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    if carbon_count < 20:
        return False, "Not a diterpenoid - requires at least 20 carbons"

    # Define SMARTS patterns for abietane core structure with variations
    core_patterns = [
        # Basic tricyclic core with methyl groups
        '[C,c]1[C,c][C,c][C,c]2[C,c]([C,c][C,c][C,c]2(C)C)[C,c]1',
        # Core with isopropyl group
        '[C,c]1[C,c][C,c]2[C,c]([C,c][C,c][C,c]2(C)C)[C,c]1C(C)C',
        # Core with oxidized positions
        '[C,c]1[C,c][C,c]2[C,c]([C,c][C,c][C,c]2(C)C(=O))[C,c]1',
        '[C,c]1[C,c][C,c]2[C,c]([C,c][C,c][C,c]2(C)CO)[C,c]1',
        # Aromatic variants
        'c1cc2[C,c][C,c][C,c]c2cc1',
        # Quinone-type structures
        'C1=CC2=C(C(=O)C1=O)[C,c][C,c][C,c]2'
    ]

    # Check for characteristic substituent patterns
    substituent_patterns = [
        # Isopropyl group
        'CC(C)[c,C]',
        # Gem-dimethyl group
        'C(C)(C)CC',
        # Common oxidized positions
        'C(=O)O[H,C]',
        'C(=O)[O-]',
        'CO[H]'
    ]

    # Check for core structure
    has_core = False
    for pattern in core_patterns:
        core = Chem.MolFromSmarts(pattern)
        if core and mol.HasSubstructMatch(core):
            has_core = True
            break

    if not has_core:
        return False, "Missing abietane core structure"

    # Check for characteristic substituents
    has_substituents = False
    for pattern in substituent_patterns:
        subst = Chem.MolFromSmarts(pattern)
        if subst and mol.HasSubstructMatch(subst):
            has_substituents = True
            break

    if has_core and has_substituents:
        # Additional check for ring connectivity
        rings = mol.GetRingInfo()
        if rings.NumRings() >= 3:  # Abietane should have at least 3 rings
            # Check for specific carbon framework
            framework = '[C,c]1[C,c][C,c][C,c]2[C,c]([C,c][C,c][C,c]2)[C,c]1'
            framework_pattern = Chem.MolFromSmarts(framework)
            if framework_pattern and mol.HasSubstructMatch(framework_pattern):
                return True, "Contains abietane core structure with characteristic substitution pattern"

    return False, "Does not match complete abietane diterpenoid structure"
# Pr=None
# Recall=0.0