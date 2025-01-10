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
    
    # Define core structural patterns for B vitamins
    patterns = {
        'B1': [  # Thiamine
            '[n+]1cscc1',  # Basic thiazolium ring
            'c1ncnc(N)c1'  # Pyrimidine part
        ],
        'B2': [  # Riboflavin
            'Cc1cc2nc3c(=O)[nH]c(=O)nc-3n(CC(O))c2cc1C'  # Isoalloxazine core
        ],
        'B3': [  # Niacin/Nicotinamide
            'c1cccnc1C(=O)[OH,N]'  # Pyridine with carboxyl or amide
        ],
        'B5': [  # Pantothenic acid
            'CC(C)(CO)C(O)C(=O)NCCC(=O)[OH]'
        ],
        'B6': [  # Pyridoxine group
            'Cc1ncc(C[OH,NH2,CHO])c(c1O)C[OH,NH2,C(=O)]'
        ],
        'B7': [  # Biotin
            'S1CC2NC(=O)NC2C1'  # Core biotin ring
        ],
        'B9': [  # Folate group
            'Nc1nc2[nH]cc(CNc3ccc(cc3)C(=O)N)nc2c(=O)[nH]1',  # Basic folate core
            'Nc1nc2NCC(CNc3ccc)nc2c(=O)[nH]1'  # Reduced forms
        ],
        'B12': [  # Cobalamin
            '[Co]',  # Must contain cobalt
            'CN1C=C2C(=C(C)C3=[N]C(=CC4=[N]C(=C(C)C5=[N]1)C)C)C' # Partial corrin pattern
        ]
    }

    # Check each vitamin pattern
    for vitamin, subpatterns in patterns.items():
        matches = 0
        for pattern in subpatterns:
            substructure = Chem.MolFromSmarts(pattern)
            if substructure is not None and mol.HasSubstructMatch(substructure):
                matches += 1
        
        # Special handling for different vitamins
        if matches > 0:
            if vitamin == 'B12':
                # Additional checks for B12
                if has_cobalamin_features(mol):
                    return True, "Contains cobalamin (B12) structure"
            elif vitamin == 'B9':
                # Check for folate characteristics
                if has_folate_features(mol):
                    return True, "Contains folate (B9) structure"
            elif matches == len(subpatterns):
                # For other vitamins, require all subpatterns to match
                return True, f"Contains vitamin {vitamin} structure"

    # Additional checks for modified forms
    if has_modified_vitamin_features(mol):
        return True, "Contains modified B vitamin structure"

    return False, "Does not match any B vitamin structural patterns"

def has_cobalamin_features(mol):
    """Check for characteristic features of B12"""
    # Must have cobalt
    if not any(atom.GetSymbol() == 'Co' for atom in mol.GetAtoms()):
        return False
    
    # Check molecular weight (cobalamins are large)
    if rdMolDescriptors.CalcExactMolWt(mol) < 1000:
        return False
    
    # Count rings (cobalamins have many)
    if rdMolDescriptors.CalcNumRings(mol) < 8:
        return False
    
    return True

def has_folate_features(mol):
    """Check for characteristic features of folates"""
    # Count key atoms
    num_n = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    num_o = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Folates typically have many nitrogens and oxygens
    if num_n < 7 or num_o < 4:
        return False
    
    # Check for characteristic p-aminobenzoyl group
    paba_pattern = Chem.MolFromSmarts('Nc1ccc(cc1)C(=O)N')
    return paba_pattern is not None and mol.HasSubstructMatch(paba_pattern)

def has_modified_vitamin_features(mol):
    """Check for common vitamin modifications"""
    # Check for phosphate groups (common in active forms)
    phosphate = Chem.MolFromSmarts('OP(=O)(O)O')
    if phosphate is not None and mol.HasSubstructMatch(phosphate):
        return True
    
    # Check for nucleotide-like features (e.g., FAD)
    nucleotide = Chem.MolFromSmarts('O1C(CO)C(O)C(O)C1n1cnc2c(N)ncnc12')
    if nucleotide is not None and mol.HasSubstructMatch(nucleotide):
        return True
    
    return False