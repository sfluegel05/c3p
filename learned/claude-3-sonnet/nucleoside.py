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
        # Sugar patterns - more flexible to catch variations
        'sugar_ring': '[CH2,CH][CR0;r5][CR0;r5][CR0;r5][OR0;r5]',
        
        # Purine patterns
        'purine_basic': 'n1cnc2c1ncnc2',
        'purine_mod': 'n1cnc2c1nc[nH]c2',
        'purine_var': '[nR1r6]1c[nR1r6]c2[nR1r6]c[nR1r6]cc12',
        
        # Pyrimidine patterns
        'pyrimidine_basic': 'n1c[nH]c(=O)[nH]c1=O',
        'pyrimidine_mod': '[nR1r6]1c(=O)[nH]c(=O)cc1',
        'pyrimidine_var': '[nR1r6]1cc[nR1r6]c(=O)[nH]1',
        
        # N-glycosidic bond - more general pattern
        'n_glycosidic': '[nR1]1[cR1][#6,#7]~[#6,#7][cR1,nR1]1[CR0]1[OR0][CR0][CR0][CR0]1',
    }
    
    # Convert patterns to RDKit molecules
    smart_patterns = {}
    for name, pattern in patterns.items():
        smart_pat = Chem.MolFromSmarts(pattern)
        if smart_pat is None:
            continue
        smart_patterns[name] = smart_pat

    # Check for sugar ring
    has_sugar = mol.HasSubstructMatch(smart_patterns['sugar_ring'])
    if not has_sugar:
        return False, "No ribose/deoxyribose/arabinose sugar found"

    # Check for nucleobase
    base_patterns = ['purine_basic', 'purine_mod', 'purine_var', 
                    'pyrimidine_basic', 'pyrimidine_mod', 'pyrimidine_var']
    found_base = False
    for base_pat in base_patterns:
        if base_pat in smart_patterns and mol.HasSubstructMatch(smart_patterns[base_pat]):
            found_base = True
            break
    
    if not found_base:
        return False, "No nucleobase pattern found"

    # Check for N-glycosidic linkage
    if not mol.HasSubstructMatch(smart_patterns['n_glycosidic']):
        # Additional check for modified N-glycosidic bonds
        ring_bonds = mol.GetSubstructMatches(smart_patterns['sugar_ring'])
        if not ring_bonds:
            return False, "No N-glycosidic bond found"
        
        # Check if any nitrogen is connected to the sugar ring
        sugar_atoms = set([atom for bond in ring_bonds for atom in bond])
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 7 and atom.IsInRing():
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIdx() in sugar_atoms:
                        found_base = True
                        break
        if not found_base:
            return False, "No valid connection between base and sugar"

    # Basic composition checks
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 7:
        return False, "Too few carbons for nucleoside"
    if n_count < 2:
        return False, "Too few nitrogens for nucleoside"
    if o_count < 2:
        return False, "Too few oxygens for nucleoside"

    # Ring checks
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:
        return False, "Nucleosides must contain at least 2 rings"

    return True, "Contains nucleobase attached to sugar via N-glycosidic bond"