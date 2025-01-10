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
    
    # Define SMARTS patterns for different B vitamins and their common forms
    patterns = {
        # B1 (Thiamine) - thiazole ring connected to pyrimidine, including charged forms
        'B1': ['[n;+]1csc([CH2]*)c1C', 'c1nc(C)nc(N)c1C[n+]1csc(CC*)c1C'],
        
        # B2 (Riboflavin) - isoalloxazine ring system
        'B2': ['Cc1cc2nc3c(=O)[nX2]c(=O)nc-3n(*)*c2cc1C'],
        
        # B3 (Niacin/Nicotinic acid) - pyridine with carboxyl
        'B3': ['O=C([OH])c1cccnc1'],
        
        # B5 (Pantothenic acid) - beta-alanine with 2,4-dihydroxy-3,3-dimethylbutanamide
        'B5': ['CC(C)(CO)[CH](O)C(=O)NCCC([OH])=O'],
        
        # B6 (Pyridoxine/Pyridoxal/Pyridoxamine) - pyridine with specific substitutions
        'B6': [
            'Cc1ncc(CO)c([CH2]*)[c]1O', # Base pattern
            'Cc1ncc(CO)c(C=O)c1O',      # Pyridoxal
            'Cc1ncc(CO)c(CN)c1O',       # Pyridoxamine
            'Cc1ncc(COP(O)(O)=O)c(*)*c1O' # Phosphate forms
        ],
        
        # B7 (Biotin) - fused rings with sulfur
        'B7': ['[H][C]12CS[C@@H](CCCCC([OH])=O)[C@@]1([H])NC(=O)N2'],
        
        # B9 (Folate) - pterin + pABA + glutamate
        'B9': [
            'Nc1nc2[nX2]cc(CNc3ccc(CC(=O)*)cc3)nc2c(=O)[nH]1',
            'Nc1nc2[nX2]c([CH2]*)c(CNc3ccc(CC(=O)*)cc3)nc2c(=O)[nH]1'
        ],
        
        # B12 (Cobalamin) - corrin ring with cobalt
        'B12': ['[Co]']
    }
    
    # Check molecular weight - B vitamins typically between 100-1500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100 or mol_wt > 2000:
        return False, "Molecular weight outside typical range for B vitamins"

    # Check each pattern
    for vitamin, subpatterns in patterns.items():
        for pattern in subpatterns:
            substructure = Chem.MolFromSmarts(pattern)
            if substructure and mol.HasSubstructMatch(substructure):
                if vitamin == 'B12':
                    # Additional check for corrin ring system for B12
                    ring_count = rdMolDescriptors.CalcNumRings(mol)
                    if ring_count >= 4 and any(atom.GetSymbol() == 'Co' for atom in mol.GetAtoms()):
                        return True, f"Matches vitamin {vitamin} (Cobalamin) structure"
                else:
                    # Verify basic composition for other vitamins
                    if check_composition(mol, vitamin):
                        return True, f"Matches vitamin {vitamin} structure"

    return False, "Does not match any B vitamin structural patterns"

def check_composition(mol, vitamin_class):
    """Helper function to verify molecular composition matches vitamin class"""
    # Count key atoms
    num_c = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    num_n = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    num_o = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    num_s = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    
    # Basic composition checks for each vitamin class
    if vitamin_class == 'B1':
        return num_c >= 12 and num_n >= 3 and num_s == 1
    elif vitamin_class == 'B2':
        return num_c >= 17 and num_n >= 4 and num_o >= 6
    elif vitamin_class == 'B3':
        return num_c >= 6 and num_n >= 1 and num_o >= 2
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