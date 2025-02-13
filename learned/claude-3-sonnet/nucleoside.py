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

    # Define and validate SMARTS patterns
    patterns = {
        # Furanose (more general pattern)
        'furanose': '[CH2,CH][CH1,CH2][CH1,CH2][CH1,CH2]O1',
        
        # Purine patterns (including modified versions)
        'purine_core': '[nr6]1[cr6][nr6][cr6]2[cr6]1[nr6][cr6][nr6][cr6]2',
        'purine_alt': '[nR1r6]1c[nR1r6]c2c1[nR1r6]c[nR1r6]c2',
        
        # Pyrimidine patterns
        'pyrimidine_core': '[nr6]1[cr6][nr6][cr6][cr6][cr6]1',
        'pyrimidine_with_O': '[nR1r6]1c[nR1r6]cc(=O)[nR1r6]c1=O',
        
        # N-glycosidic bond patterns (more general)
        'n_glycosidic': '[nr6]C1[O][C@@H,C@H,CH]([CH2,CH])[CH1,CH2][CH1,CH2]1',
        'n_glycosidic_alt': '[nR1]C1OC(CO)CC1'
    }
    
    # Validate and convert patterns
    smart_patterns = {}
    for name, pattern in patterns.items():
        smart_pat = Chem.MolFromSmarts(pattern)
        if smart_pat is None:
            continue
        smart_patterns[name] = smart_pat

    # Check for furanose (sugar ring)
    if not any(mol.HasSubstructMatch(pat) for pat in 
              [smart_patterns['furanose']] if 'furanose' in smart_patterns):
        return False, "No ribose/deoxyribose sugar found"

    # Check for nucleobase
    base_patterns = ['purine_core', 'purine_alt', 'pyrimidine_core', 'pyrimidine_with_O']
    found_base = False
    for base_pat in base_patterns:
        if base_pat in smart_patterns and mol.HasSubstructMatch(smart_patterns[base_pat]):
            found_base = True
            break
    
    if not found_base:
        return False, "No nucleobase pattern found"

    # Check for N-glycosidic bond
    glycosidic_patterns = ['n_glycosidic', 'n_glycosidic_alt']
    found_glycosidic = False
    for gly_pat in glycosidic_patterns:
        if gly_pat in smart_patterns and mol.HasSubstructMatch(smart_patterns[gly_pat]):
            found_glycosidic = True
            break
            
    if not found_glycosidic:
        return False, "No N-glycosidic bond found"

    # Basic size and composition checks
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight {mol_wt:.1f} too small for nucleoside"
    
    # Count atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 8:
        return False, "Too few carbons for nucleoside"
    if n_count < 2:
        return False, "Too few nitrogens for nucleoside"
    if o_count < 3:
        return False, "Too few oxygens for nucleoside"

    # Check ring count
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:
        return False, "Too few rings for nucleoside structure"
    
    # Additional check for aromatic atoms (nucleobases are typically aromatic)
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    if aromatic_atoms < 4:
        return False, "Too few aromatic atoms for typical nucleoside"

    return True, "Contains nucleobase attached to ribose/deoxyribose via N-glycosidic bond"