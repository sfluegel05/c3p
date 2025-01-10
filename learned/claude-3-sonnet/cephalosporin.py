"""
Classifies: CHEBI:23066 cephalosporin
"""
"""
Classifies: cephalosporin antibiotics
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_cephalosporin(smiles: str):
    """
    Determines if a molecule is a cephalosporin based on its SMILES string.
    Cephalosporins contain a beta-lactam ring fused to a 6-membered dihydrothiazine ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cephalosporin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Split into fragments (for salts/hydrates) and get largest fragment
    fragments = Chem.GetMolFrags(mol, asMols=True)
    if fragments:
        mol = max(fragments, key=lambda m: m.GetNumAtoms())

    # Core cephalosporin patterns - multiple SMARTS to catch all variations
    core_patterns = [
        # Basic cephalosporin core (more general)
        'S1CC2N(C(=O)C2*)C(=C1*)*',
        # Alternative with explicit connection points
        'S1CC2N(C(=O)C2NC*)C(=C1C*)C(=O)[O,N]',
        # More generic pattern allowing for variations
        'S1[CH2,CH1]C2N([*])C(=O)C2[*]C1[*]',
        # Pattern for saturated versions
        'S1[CH2][CH1]2N([*])C(=O)[CH1]2[*]C1[*]'
    ]
    
    found_core = False
    matching_pattern = None
    for pattern in core_patterns:
        core = Chem.MolFromSmarts(pattern)
        if core and mol.HasSubstructMatch(core):
            found_core = True
            matching_pattern = pattern
            break
            
    if not found_core:
        return False, "Missing cephalosporin core structure"

    # Check for carboxylic acid/carboxylate group
    carboxyl_patterns = [
        'C(=O)[OH]',      # carboxylic acid
        'C(=O)[O-]',      # carboxylate
        'C(=O)O[Na,K]',   # metal salts
        'C(=O)O'          # generic carboxyl
    ]
    
    has_carboxyl = False
    for pattern in carboxyl_patterns:
        p = Chem.MolFromSmarts(pattern)
        if p and mol.HasSubstructMatch(p):
            has_carboxyl = True
            break
            
    if not has_carboxyl:
        return False, "Missing carboxylic acid group"

    # Check for characteristic substituents commonly found in cephalosporins
    substituent_patterns = {
        'aminothiazole': 'c1csc(N)n1',
        'tetrazole': 'n1nnn[nH,c]1',
        'oxime': '[CX3]=[NX2][OX2][*]',
        'acetoxymethyl': 'CC(=O)OC',
        'amino': '[NH2]',
        'amide': 'NC(=O)',
        'thiadiazole': 'c1nncs1'
    }
    
    found_substituents = []
    for name, pattern in substituent_patterns.items():
        p = Chem.MolFromSmarts(pattern)
        if p and mol.HasSubstructMatch(p):
            found_substituents.append(name)

    # Additional structural checks
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    
    # Count S atoms (should have at least one in core)
    s_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16])
    
    # Count N atoms (should have at least 2 - one in beta-lactam, one in side chain)
    n_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7])

    if num_rings < 2 or s_count < 1 or n_count < 2:
        return False, "Missing essential structural elements for cephalosporin"

    # Build reason string
    reason = "Identified as cephalosporin:\n"
    reason += "- Contains beta-lactam fused to dihydrothiazine ring system\n"
    if has_carboxyl:
        reason += "- Contains carboxylic acid/carboxylate group\n"
    if found_substituents:
        reason += f"- Found characteristic substituents: {', '.join(found_substituents)}\n"
    reason += f"- Structure contains {num_rings} rings, {s_count} sulfur atoms, {n_count} nitrogen atoms"

    return True, reason