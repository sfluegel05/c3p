"""
Classifies: CHEBI:75769 B vitamin
"""
"""
Classifies: B vitamin compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_B_vitamin(smiles: str):
    """
    Determines if a molecule is a B vitamin based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a B vitamin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Exclude pattern - common non-vitamin structures that might give false positives
    exclude_patterns = [
        'F[c]1[c][n]2[c](=O)[c](C(=O)O)[c][n]', # Fluoroquinolone core
        '[$(C1=CC=C2N(C=C(C(=O)O)C(=O)C2=C1)C)]' # Quinolone core
    ]
    
    for pattern in exclude_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return False, "Contains structural features typical of non-vitamin compounds"

    # Define SMARTS patterns for different B vitamins and their common forms
    patterns = {
        # B1 (Thiamine) - thiazole ring connected to pyrimidine
        'B1': [
            '[n;+]1csc(CC*)c1C',  # Thiazolium core
            'c1nc(C)nc(N)c1C[n+]1csc(CC*)c1C', # Complete thiamine
            'c1nc(C)nc(N)c1C[n+]1csc(CCOP*)c1C' # Phosphorylated forms
        ],
        
        # B2 (Riboflavin) - isoalloxazine ring system and derivatives
        'B2': [
            'Cc1cc2nc3c(=O)[nX2]c(=O)nc-3n(CC(O)C(O)C(O)*)c2cc1C', # Riboflavin
            'Cc1cc2nc3c(=O)[n-]c(=O)nc-3n(CC(O)C(O)C(O)*)c2cc1C',  # Reduced form
            'Cc1cc2nc3c(=O)[nX2]c(=O)nc-3n(CC(O)C(O)C(O)COP*)c2cc1C' # FMN/FAD core
        ],
        
        # B3 (Niacin) - pyridine with specific substitution pattern
        'B3': [
            '[$(O=C(O)c1cccnc1):1][$([H,O-]):2]', # Nicotinic acid
            '[$(O=C(O)c1ccc[n;H]c1):1][$([H,O-]):2]', # Reduced form
            'CN1C=CC(=CC1=O)C(=O)[O;H,-]' # NAD-related core
        ],
        
        # B5 (Pantothenic acid) - specific chain with stereochemistry
        'B5': [
            'CC(C)(CO)[C@@H](O)C(=O)NCCC([OH,O-])=O',
            'CC(C)(COP*)C(O)C(=O)NCCC([OH,O-])=O' # Phosphorylated form
        ],
        
        # B6 group - pyridoxine/pyridoxal/pyridoxamine
        'B6': [
            'Cc1ncc(CO)c(C[NH2,NH3+])c1O', # Pyridoxamine
            'Cc1ncc(CO)c(C=O)c1O',         # Pyridoxal
            'Cc1ncc(CO)c(CO)c1O',          # Pyridoxine
            'Cc1ncc(COP(O)(O)=O)c(C*)c1O', # Phosphate forms
            'Cc1ncc(CO)c(C(=O)[O;H,-])c1O' # Acid form
        ],
        
        # B7 (Biotin)
        'B7': [
            '[H][C@]12CS[C@@H](CCCCC([OH,O-])=O)[C@@]1([H])NC(=O)N2',
            '[H][C@]12CS[C@@H](CCCCCC([OH,O-])=O)[C@@]1([H])NC(=O)N2' # Homobiotin
        ],
        
        # B9 (Folate) group - including reduced forms
        'B9': [
            'Nc1nc2N[CH2,CH]C(CNc3ccc(CC(=O)N[CH]CC([OH,O-])=O)cc3)Nc2c(=O)[nH]1', # THF
            'Nc1nc2N=CC(CNc3ccc(CC(=O)N[CH]CC([OH,O-])=O)cc3)Nc2c(=O)[nH]1',       # DHF
            'Nc1nc2ncc(CNc3ccc(CC(=O)N[CH]CC([OH,O-])=O)cc3)nc2c(=O)[nH]1',        # Folate
            '[CH3,CH2OH,CHO,CH=NH]N1[CH]2CNc3nc(N)[nH]c(=O)c3N(*)C2CN1*'           # Modified forms
        ],
        
        # B12 (Cobalamin) - corrin ring with cobalt
        'B12': [
            '[Co]',  # Must contain cobalt
            '[Co]N4' # Cobalt with corrin coordination
        ]
    }
    
    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100 or mol_wt > 2000:
        return False, "Molecular weight outside typical range for B vitamins"

    # Check each pattern
    for vitamin, subpatterns in patterns.items():
        for pattern in subpatterns:
            substructure = Chem.MolFromSmarts(pattern)
            if substructure and mol.HasSubstructMatch(substructure):
                if vitamin == 'B12':
                    # Additional checks for B12
                    if has_corrin_system(mol):
                        return True, f"Matches vitamin {vitamin} (Cobalamin) structure"
                else:
                    # Verify composition and additional features
                    if check_composition(mol, vitamin):
                        return True, f"Matches vitamin {vitamin} structure"

    return False, "Does not match any B vitamin structural patterns"

def has_corrin_system(mol):
    """Check for characteristic corrin ring system of B12"""
    # Count rings and check for characteristic coordination
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    has_co = any(atom.GetSymbol() == 'Co' for atom in mol.GetAtoms())
    # Look for characteristic nitrogen coordination around cobalt
    co_coordination = Chem.MolFromSmarts('[Co]~N~C~C~N')
    
    return ring_count >= 4 and has_co and mol.HasSubstructMatch(co_coordination)

def check_composition(mol, vitamin_class):
    """Helper function to verify molecular composition matches vitamin class"""
    # Count key atoms
    num_c = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    num_n = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    num_o = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    num_s = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    num_p = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    # Adjusted composition checks
    if vitamin_class == 'B1':
        return num_c >= 12 and num_n >= 3 and num_s == 1
    elif vitamin_class == 'B2':
        return num_c >= 15 and num_n >= 4 and num_o >= 4
    elif vitamin_class == 'B3':
        return 6 <= num_c <= 12 and num_n >= 1 and num_o >= 2 and num_s == 0
    elif vitamin_class == 'B5':
        return num_c >= 9 and num_n >= 1 and num_o >= 4
    elif vitamin_class == 'B6':
        return num_c >= 8 and num_n >= 1 and num_o >= 2
    elif vitamin_class == 'B7':
        return num_c >= 10 and num_n >= 2 and num_o >= 2 and num_s == 1
    elif vitamin_class == 'B9':
        return num_c >= 19 and num_n >= 7 and num_o >= 6
    
    return True

__metadata__ = {
    'chemical_class': {
        'name': 'B vitamin',
        'definition': 'Any member of the group of eight water-soluble vitamins '
                     'originally thought to be a single compound (vitamin B) that '
                     'play important roles in cell metabolism.',
        'subclasses': ['B1', 'B2', 'B3', 'B5', 'B6', 'B7', 'B9', 'B12']
    }
}