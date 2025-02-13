"""
Classifies: CHEBI:17334 penicillin
"""
"""
Classifies: CHEBI:17334 penicillin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_penicillin(smiles: str):
    """
    Determines if a molecule is a penicillin based on its SMILES string.
    Penicillins have a characteristic penam core (beta-lactam fused to thiazolidine ring)
    with specific substituents: two methyls at position 2, carboxylate at position 3,
    and carboxamido at position 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a penicillin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for penam core structure (beta-lactam fused to thiazolidine ring)
    # Pattern includes the key stereocenters and ring fusion
    penam_patterns = [
        # Pattern 1: More specific with stereochemistry
        '[H][C@]12S[C@@](C)(C)[C@@H](N1)C(=O)N2',
        # Pattern 2: Alternative representation
        'S1[C@]2([H])[C@@H](N[C@@H]1C(C)(C))[C(=O)N2',
        # Pattern 3: More general pattern
        'S1C2NC(=O)C2N1'
    ]
    
    core_found = False
    for pattern in penam_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            core_found = True
            break
            
    if not core_found:
        return False, "No penam core structure found"

    # Check for carboxylate group at position 3 (more flexible pattern)
    carboxylate_patterns = [
        'S[C@](C)(C)[C@@H](N)C(=O)[O,OH]',  # Free acid form
        'S[C@](C)(C)[C@@H](N)C(=O)O[C,H]'   # Ester form
    ]
    
    carboxylate_found = False
    for pattern in carboxylate_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            carboxylate_found = True
            break
            
    if not carboxylate_found:
        return False, "Missing carboxylate group at position 3"

    # Check for two methyl groups at position 2
    dimethyl_pattern = Chem.MolFromSmarts('SC(C)(C)[C,H]')
    if not mol.HasSubstructMatch(dimethyl_pattern):
        return False, "Missing two methyl groups at position 2"

    # Check for carboxamido group at position 6
    carboxamido_pattern = Chem.MolFromSmarts('[NH][CH]1C(=O)N2')
    if not mol.HasSubstructMatch(carboxamido_pattern):
        return False, "Missing carboxamido group at position 6"

    # Verify overall composition
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    if s_count != 1:
        return False, f"Incorrect number of sulfur atoms (found {s_count}, expected 1)"

    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 2:
        return False, f"Too few nitrogen atoms (found {n_count}, expected at least 2)"

    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:
        return False, f"Too few oxygen atoms (found {o_count}, expected at least 3)"

    # Check ring system
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:
        return False, "Missing required bicyclic system"

    return True, "Contains penam core with correct substituents"