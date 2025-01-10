"""
Classifies: CHEBI:63534 monoamine
"""
"""
Classifies: monoamine compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.
    A monoamine has one amino group connected to an aromatic ring by a two-carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES and handle salt forms
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get the largest fragment (in case of salts)
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) > 1:
        mol = max(frags, key=lambda m: m.GetNumAtoms())

    # Basic structure requirements
    if mol.GetNumAtoms() < 8:  # Minimum size for a monoamine
        return False, "Molecule too small for monoamine structure"
    
    if mol.GetNumAtoms() > 50:  # Maximum reasonable size
        return False, "Molecule too large for typical monoamine"

    # Check for aromatic ring
    aromatic_ring = Chem.MolFromSmarts("a1aaaaa1")  # 6-membered aromatic ring
    if not mol.HasSubstructMatch(aromatic_ring):
        return False, "No aromatic ring found"

    # Define core monoamine patterns more precisely
    monoamine_patterns = [
        # Basic ethylamine patterns with explicit connection to aromatic ring
        "a[CH2][CH2][NH2]",  # Simple primary amine
        "a[CH2][CH2][NH][CH3]",  # Secondary amine (methyl)
        "a[CH2][CH2][NH]C",  # Secondary amine (any alkyl)
        "a[CH2][CH2][N]([CH3])[CH3]",  # Tertiary amine (dimethyl)
        
        # Beta-hydroxyl patterns
        "a[CH2][CH]([OH])[NH2]",
        "a[CH2][CH]([OH])[NH][CH3]",
        "a[CH2][CH]([OH])[N]([CH3])[CH3]",
        
        # Charged variants
        "a[CH2][CH2][NH3+]",
        "a[CH2][CH2][NH2+][CH3]",
        "a[CH2][CH2][NH+]([CH3])[CH3]",
        
        # Alpha-hydroxyl patterns
        "a[CH]([OH])[CH2][NH2]",
        "a[CH]([OH])[CH2][NH][CH3]",
        "a[CH]([OH])[CH2][N]([CH3])[CH3]"
    ]

    found_pattern = False
    for pattern in monoamine_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt and mol.HasSubstructMatch(patt):
            found_pattern = True
            break
            
    if not found_pattern:
        return False, "No valid monoamine pattern found"

    # Count total amine groups (including charged)
    amine_patterns = [
        "[NH2]", "[NH][CH3]", "[N]([CH3])[CH3]",  # Neutral amines
        "[NH3+]", "[NH2+][CH3]", "[NH+]([CH3])[CH3]"  # Charged amines
    ]
    
    total_amines = 0
    for pattern in amine_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt:
            total_amines += len(mol.GetSubstructMatches(patt))

    if total_amines == 0:
        return False, "No amine groups found"
    if total_amines > 2:  # Allow maximum two amine groups
        return False, "Too many amine groups"

    # Additional checks for common monoamine features
    aromatic_oh = Chem.MolFromSmarts("aO[H]")
    has_aromatic_oh = mol.HasSubstructMatch(aromatic_oh)
    
    # Check for problematic features that might indicate non-monoamine
    problematic_features = [
        (Chem.MolFromSmarts("[N]1[C]=[N][C]=[N][C]1"), "Contains tetrazole"),
        (Chem.MolFromSmarts("[N]1[C]=[O][C]="), "Contains oxazole"),
        (Chem.MolFromSmarts("[N]1[N]=[N][N]="), "Contains tetrazine"),
        (Chem.MolFromSmarts("C(=O)O[H]"), "Contains free carboxylic acid")
    ]
    
    for pattern, reason in problematic_features:
        if pattern and mol.HasSubstructMatches(pattern):
            matches = mol.GetSubstructMatches(pattern)
            if len(matches) > 1:  # Allow one instance
                return False, f"Multiple instances of {reason}"

    detail = "Contains aromatic ring with"
    if has_aromatic_oh:
        detail += " hydroxyl group(s) and"
    detail += " two-carbon chain connected to amine group"
    
    return True, detail