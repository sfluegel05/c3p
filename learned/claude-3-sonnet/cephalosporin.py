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
    
    # Parse SMILES - consider largest fragment for salts
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Split into fragments (for salts/hydrates) and get largest fragment
    fragments = Chem.GetMolFrags(mol, asMols=True)
    if fragments:
        mol = max(fragments, key=lambda m: m.GetNumAtoms())

    # Core cephalosporin patterns - trying multiple SMARTS to catch variations
    core_patterns = [
        # Basic ceph core with explicit stereochemistry
        '[H][C@@]12SCC=C(N1C(=O)[C@H]2[*])[*]',
        # Alternative core pattern
        'S1CC=C2N([*])C(=O)[C@H]2N1',
        # More generic pattern
        'S1CC=C2N(C(=O)C2[*])C1=C[*]'
    ]
    
    found_core = False
    for pattern in core_patterns:
        core = Chem.MolFromSmarts(pattern)
        if core and mol.HasSubstructMatch(core):
            found_core = True
            break
            
    if not found_core:
        return False, "Missing cephalosporin core structure"

    # Check for carboxylic acid/carboxylate group
    carboxyl_patterns = [
        'C(=O)[OH]',      # carboxylic acid
        'C(=O)[O-]',      # carboxylate
        'C(=O)O[Na]',     # sodium salt
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

    # Check for characteristic substituents
    substituent_patterns = {
        'aminothiazole': 'c1csc(N)n1',
        'tetrazole': 'n1nnn[nH]1',
        'methyltetrazole': 'Cn1nnnc1',
        'oxime': 'C=NOC',
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

    # Verify basic properties
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 15:
        return False, "Molecule too small to be a cephalosporin"

    # Count rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:
        return False, "Insufficient ring count for cephalosporin structure"

    # Build reason string
    reason = "Identified as cephalosporin due to:\n"
    reason += "- Contains fused beta-lactam-dihydrothiazine ring system\n"
    reason += "- Contains carboxylic acid/carboxylate group\n"
    if found_substituents:
        reason += f"- Found characteristic substituents: {', '.join(found_substituents)}\n"
    reason += f"- Contains appropriate molecular complexity ({num_atoms} atoms, {ring_info.NumRings()} rings)"

    return True, reason