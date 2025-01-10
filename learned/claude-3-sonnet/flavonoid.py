"""
Classifies: CHEBI:47916 flavonoid
"""
"""
Classifies: CHEBI:47916 flavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    A flavonoid has a 1-benzopyran core with an aryl substituent at position 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic size and composition check
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 15:
        return False, "Molecule too small to be a flavonoid"

    # Core patterns for different flavonoid classes
    core_patterns = [
        # Basic flavone/flavonol core (more flexible)
        "[#6]1=C([#6])Oc2c([#6,#1])c([#6,#1,#8])c([#6,#1,#8])c([#6,#1,#8])c2[#6]1=O",
        
        # Flavanone core
        "[#6]1[#6]Oc2c([#6,#1])c([#6,#1,#8])c([#6,#1,#8])c([#6,#1,#8])c2[#6]1=O",
        
        # Isoflavone core
        "O=C1C(=C)Oc2c([#6,#1])c([#6,#1,#8])c([#6,#1,#8])c([#6,#1,#8])c21",
        
        # Anthocyanidin core (includes charged species)
        "[#6]1=[#6+]([#6])Oc2c([#6,#1])c([#6,#1,#8])c([#6,#1,#8])c([#6,#1,#8])c21",
        
        # More general chromane/chromene core
        "[#6]1[#6]Oc2c([#6,#1])c([#6,#1,#8])c([#6,#1,#8])c([#6,#1,#8])c21",
        
        # Flavan core
        "[#6]1[#6][#6]Oc2c([#6,#1])c([#6,#1,#8])c([#6,#1,#8])c([#6,#1,#8])c21"
    ]

    has_core = False
    for pattern in core_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat is not None and mol.HasSubstructMatch(pat):
            has_core = True
            break

    if not has_core:
        return False, "No flavonoid core structure found"

    # Check for aromatic/phenyl substituent
    aryl_patterns = [
        # General aryl group patterns
        "c1ccccc1",
        "c1cc(O)ccc1",
        "c1cc(O)cc(O)c1",
        
        # More specific patterns for common substitution
        "c1c(O)c(O)ccc1",
        "c1c(OC)c(O)ccc1"
    ]

    has_aryl = False
    for pattern in aryl_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat is not None and mol.HasSubstructMatch(pat):
            has_aryl = True
            break

    if not has_aryl:
        return False, "Missing required aryl substituent"

    # Check for typical flavonoid characteristics
    characteristic_patterns = [
        # Hydroxyl groups
        "[OX2H1]",
        # Methoxy groups
        "[OX2][CH3]",
        # Sugar linkages (more flexible)
        "[OX2][CH]([OH])[CH]([OH])[CH]",
        # Carbonyl groups
        "[CX3](=[OX1])",
        # Glycosidic linkages
        "[OX2]([CH]1[OH0,OH1])[CH]([OH0,OH1])[CH]([OH0,OH1])[CH]([OH0,OH1])[CH]([OH0,OH1])[CH]1"
    ]

    characteristic_count = 0
    for pattern in characteristic_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat is not None:
            characteristic_count += len(mol.GetSubstructMatches(pat))

    if characteristic_count < 2:
        return False, "Insufficient number of characteristic flavonoid substituents"

    # Ring count check
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:
        return False, "Insufficient number of rings"

    # Additional check for fused ring systems
    if mol.GetNumBonds() < 20:  # Most flavonoids have at least this many bonds
        return False, "Structure too simple for a flavonoid"

    return True, "Contains flavonoid core structure with appropriate substituents"