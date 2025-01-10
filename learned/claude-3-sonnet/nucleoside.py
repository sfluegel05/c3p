"""
Classifies: CHEBI:33838 nucleoside
"""
"""
Classifies: CHEBI:33838 nucleoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside(smiles: str):
    """
    Determines if a molecule is a nucleoside based on its SMILES string.
    A nucleoside consists of a nucleobase attached to ribose/deoxyribose via N-glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns
    patterns = {
        # Furanose sugar patterns - more flexible to catch variations
        'sugar_ring': '[#6]1[#6][#6][#6][#8]1', # Basic furanose
        'ribose': '[#6]1([#6][#8])[#8][#6]([#6]([#8])[#6]1[#8])',  # More specific ribose pattern
        
        # Nucleobase patterns
        'purine_core': 'c1ncnc2[nX2]cnc12',  # Basic purine scaffold
        'purine_mod': 'c1nc([*,#7,#8])nc2[nX2]cnc12',  # Modified purine
        'pyrimidine_core': 'c1cn([*])[cX3][cX3,nX2]n1',  # Basic pyrimidine
        'pyrimidine_mod': 'c1c([*,#7,#8])n([*])[cX3][cX3,nX2]n1',  # Modified pyrimidine
        
        # N-glycosidic bond - simplified
        'n_glycosidic': '[#6]1[#8][#6][#6][#6]1[NX3]',  # Sugar-N connection
    }
    
    # Convert patterns to RDKit molecules
    smart_patterns = {}
    for name, pattern in patterns.items():
        smart_pat = Chem.MolFromSmarts(pattern)
        if smart_pat is None:
            continue
        smart_patterns[name] = smart_pat

    # Check for sugar ring (either basic furanose or specific ribose pattern)
    has_sugar = False
    if mol.HasSubstructMatch(smart_patterns['sugar_ring']):
        has_sugar = True
    if mol.HasSubstructMatch(smart_patterns['ribose']):
        has_sugar = True
    
    if not has_sugar:
        return False, "No ribose/deoxyribose sugar found"

    # Check for nucleobase
    has_base = False
    base_patterns = ['purine_core', 'purine_mod', 'pyrimidine_core', 'pyrimidine_mod']
    for pattern in base_patterns:
        if pattern in smart_patterns and mol.HasSubstructMatch(smart_patterns[pattern]):
            has_base = True
            break

    if not has_base:
        return False, "No nucleobase pattern found"

    # Check for N-glycosidic linkage
    if not mol.HasSubstructMatch(smart_patterns['n_glycosidic']):
        return False, "No N-glycosidic bond found"

    # Additional structural checks
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:
        return False, "Must contain at least 2 rings (sugar + base)"

    # Count key atoms
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if n_count < 1:
        return False, "Must contain at least 1 nitrogen"
    if o_count < 3:
        return False, "Must contain at least 3 oxygens"

    # Check ring sizes
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    if 5 not in ring_sizes:
        return False, "Must contain 5-membered sugar ring"
    
    # Calculate basic properties
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 200 or mw > 800:
        return False, "Molecular weight outside typical range for nucleosides"

    return True, "Contains nucleobase attached to sugar via N-glycosidic bond"