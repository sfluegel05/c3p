"""
Classifies: CHEBI:46633 carbapenems
"""
"""
Classifies: CHEBI:23066 carbapenem
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_carbapenems(smiles: str):
    """
    Determines if a molecule is a carbapenem based on its SMILES string.
    Carbapenems have a beta-lactam ring fused to a 5-membered ring with
    various substitutions at positions 2, 3, and 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbapenem, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Essential core patterns - beta-lactam fused to 5-membered ring
    core_patterns = [
        # Basic beta-lactam fusion pattern (more flexible)
        Chem.MolFromSmarts("[#6]1~[#6]~[#6]~[#7]2~[#6]~[#6]12"),
        # Alternative representation
        Chem.MolFromSmarts("[#6]2~[#6]~[#6]~[#7]1~[#6](=[O])~[#6]12"),
        # Pattern with carbonyl
        Chem.MolFromSmarts("[#6]1~[#6]~[#6]~[#7]2~[#6](=O)~[#6]12")
    ]
    
    core_match = False
    for pattern in core_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            core_match = True
            break
            
    if not core_match:
        return False, "Missing carbapenem core structure (fused beta-lactam ring system)"

    # Check for ring sizes
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    if not (4 in ring_sizes and 5 in ring_sizes):
        return False, "Must contain both 4-membered (beta-lactam) and 5-membered rings"

    # Look for characteristic groups
    patterns = {
        'carboxyl': [
            Chem.MolFromSmarts("C(=O)O"),
            Chem.MolFromSmarts("C(=O)[O-]")
        ],
        'sulfur': [
            Chem.MolFromSmarts("[#6]~[#16]"),  # C-S bond
            Chem.MolFromSmarts("[#6]S[#6]")    # C-S-C
        ],
        'hydroxy': [
            Chem.MolFromSmarts("CO"),          # Any hydroxyl
            Chem.MolFromSmarts("CC(O)C"),      # Secondary alcohol
            Chem.MolFromSmarts("CC(O)")        # Terminal hydroxyl
        ],
        'beta_lactam': [
            Chem.MolFromSmarts("[#7]1[#6](=O)[#6][#6]1"),  # Beta-lactam ring
            Chem.MolFromSmarts("[#7]1[#6](=O)[#6][#6]1")   # Alternative representation
        ]
    }

    # Check for presence of key features
    features = []
    for feature, pattern_list in patterns.items():
        for pattern in pattern_list:
            if pattern is not None and mol.HasSubstructMatch(pattern):
                features.append(feature)
                break

    # Must have beta-lactam and at least one other feature
    if 'beta_lactam' not in features:
        return False, "Missing beta-lactam ring"
    
    if len(set(features)) < 2:  # Using set to count unique features
        return False, "Insufficient characteristic carbapenem features"

    # Additional checks
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    mol_wt = Chem.Descriptors.ExactMolWt(mol)

    if n_count < 1:
        return False, "Must contain at least one nitrogen atom"
    if o_count < 2:
        return False, "Must contain at least two oxygen atoms"
    if mol_wt < 150:  # Lowered threshold to catch simpler core structures
        return False, "Molecular weight too low for carbapenem"

    # Success case
    unique_features = set(features) - {'beta_lactam'}  # Remove beta_lactam from feature list
    reason = f"Contains carbapenem core structure with {', '.join(unique_features)}"
    return True, reason