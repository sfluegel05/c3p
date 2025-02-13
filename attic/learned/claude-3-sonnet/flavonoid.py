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

    # Core patterns for different flavonoid classes - simplified and more general
    core_patterns = [
        # Basic benzopyran core (most general)
        "O1c2ccccc2CC1",
        
        # Flavone/flavonol core (more general)
        "O1c2ccccc2C(=O)C1",
        
        # Isoflavone core
        "O1c2ccccc2C(=O)C=C1",
        
        # Flavanone core
        "O1CCc2ccccc2C1=O",
        
        # Anthocyanidin core (includes charged form)
        "O1c2ccccc2C=[O+]C1",
        
        # Chalcone-type pattern
        "O=CC(=O)c1ccccc1O",
        
        # More specific but common patterns
        "O1c2c(O)cc(O)cc2OC(c2ccccc2)C1=O",
        "O1c2c(O)cc(O)cc2OC(c2ccccc2)C1",
    ]

    has_core = False
    matched_pattern = None
    for pattern in core_patterns:
        pat = Chem.MolFromSmiles(pattern)
        if pat is not None and mol.HasSubstructMatch(pat):
            has_core = True
            matched_pattern = pattern
            break

    if not has_core:
        # Try more general SMARTS patterns if SMILES patterns fail
        smarts_patterns = [
            # Very general benzopyran core
            "O1[#6]~2~[#6]~[#6]~[#6]~[#6]~[#6]2[#6][#6]1",
            # General chromone core
            "O1[#6]~2~[#6]~[#6]~[#6]~[#6]~[#6]2[#6](=O)[#6]1",
        ]
        for pattern in smarts_patterns:
            pat = Chem.MolFromSmarts(pattern)
            if pat is not None and mol.HasSubstructMatch(pat):
                has_core = True
                matched_pattern = pattern
                break

    if not has_core:
        return False, "No flavonoid core structure found"

    # Check for aromatic/phenyl substituent - more general patterns
    aryl_patterns = [
        "c1ccccc1",  # Basic phenyl
        "c1cccc([#8,#7,#6])c1",  # Substituted phenyl
        "c1cc([#8,#7])cc([#8,#7])c1",  # Di-substituted phenyl
    ]

    has_aryl = False
    for pattern in aryl_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat is not None and mol.HasSubstructMatch(pat):
            has_aryl = True
            break

    # Characteristic flavonoid features
    characteristic_patterns = [
        "[OX2H1]",  # Hydroxyl
        "[OX2][CH3]",  # Methoxy
        "C(=O)",  # Carbonyl
        "O[CH]1O[CH][CH][CH][CH][CH]1",  # Sugar ring
        "O[CH]1O[CH][CH][CH][CH]1",  # Furanose sugar
    ]

    characteristic_count = 0
    for pattern in characteristic_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat is not None:
            matches = len(mol.GetSubstructMatches(pat))
            characteristic_count += matches

    # Ring count check
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:
        return False, "Insufficient number of rings"

    # Classification logic
    if has_aryl and characteristic_count >= 1:
        return True, "Contains flavonoid core structure with appropriate substituents"
    elif has_core and characteristic_count >= 2:
        return True, "Contains flavonoid core with characteristic substitution pattern"
    else:
        return False, "Missing required structural features for flavonoid classification"