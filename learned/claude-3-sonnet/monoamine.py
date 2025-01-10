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
    
    # Check for aromatic ring
    aromatic_pattern = Chem.MolFromSmarts("a1aaaaa1")  # 6-membered aromatic ring
    if not mol.HasSubstructMatch(aromatic_pattern):
        return False, "No aromatic ring found"

    # Define monoamine patterns more precisely
    # Key features:
    # - Must have exactly two carbons between aromatic ring and amine
    # - Allow for substitutions on the carbons
    # - Handle both charged and uncharged forms
    # - Exclude amides and ring nitrogens
    monoamine_patterns = [
        # Basic ethylamine pattern
        "[aR1]!@[CH2][CH2][NX3;!R;!$(NC=O)]",
        # Beta-hydroxyl pattern (like in adrenaline)
        "[aR1]!@[CH2][CH1]([OH1])[NX3;!R;!$(NC=O)]",
        # Alpha-hydroxyl pattern
        "[aR1]!@[CH1]([OH1])[CH2][NX3;!R;!$(NC=O)]",
        # Methylated variants
        "[aR1]!@[CH2][CH2][NX3H2;!R;!$(NC=O)]C",
        "[aR1]!@[CH2][CH2][NX3H1;!R;!$(NC=O)](C)C",
        # Charged variants
        "[aR1]!@[CH2][CH2][NX4H3+;!R]",
        "[aR1]!@[CH2][CH2][NX4H2+;!R]C",
        # Additional substituted patterns
        "[aR1]!@[CH2][CH1]([*])[NX3;!R;!$(NC=O)]",
        "[aR1]!@[CH1]([*])[CH2][NX3;!R;!$(NC=O)]"
    ]

    found_valid_pattern = False
    for pattern in monoamine_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt and mol.HasSubstructMatch(patt):
            found_valid_pattern = True
            break
            
    if not found_valid_pattern:
        return False, "No valid two-carbon chain connecting amine to aromatic ring"

    # Additional validations
    
    # Count non-aromatic ring nitrogens (exclude molecules with too many amines)
    amine_pattern = Chem.MolFromSmarts("[NX3;!R;!$(NC=O)]")
    charged_amine_pattern = Chem.MolFromSmarts("[NX4H+;!R]")
    
    amine_count = len(mol.GetSubstructMatches(amine_pattern))
    charged_amine_count = len(mol.GetSubstructMatches(charged_amine_pattern))
    total_amines = amine_count + charged_amine_count
    
    if total_amines > 3:  # allowing up to 3 for some derivatives
        return False, "Too many amine groups"

    # Exclude if nitrogen is part of a ring system
    ring_n_pattern = Chem.MolFromSmarts("[NR]")
    if mol.HasSubstructMatch(ring_n_pattern):
        # Check if ALL nitrogens are ring nitrogens
        n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
        ring_n_count = len(mol.GetSubstructMatches(ring_n_pattern))
        if ring_n_count == n_count:
            return False, "All nitrogens are part of rings"

    # Look for common monoamine features (hydroxyl groups on aromatic ring)
    aromatic_oh_pattern = Chem.MolFromSmarts("aO[H]")
    has_aromatic_oh = mol.HasSubstructMatch(aromatic_oh_pattern)
    
    detail = "monoamine with"
    if has_aromatic_oh:
        detail += " hydroxylated"
    detail += " aromatic ring and two-carbon amine linker"
    
    return True, f"Confirmed {detail}"