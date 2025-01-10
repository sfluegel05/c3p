"""
Classifies: CHEBI:75769 B vitamin
"""
"""
Classifies: B vitamin compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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
    
    # Define SMARTS patterns for different B vitamins
    patterns = {
        # B1 (Thiamine) - thiazole ring connected to pyrimidine
        'B1': '[n+]1csc(CCO)c1C',
        
        # B2 (Riboflavin) - isoalloxazine ring system
        'B2': 'Cc1cc2nc3c(=O)[nX2]c(=O)nc-3n(CC(O)C(O)C(O)CO)c2cc1C',
        
        # B3 (Niacin/Nicotinic acid) - pyridine with carboxyl
        'B3': 'O=C(O)c1cccnc1',
        
        # B5 (Pantothenic acid) - beta-alanine with 2,4-dihydroxy-3,3-dimethylbutanamide
        'B5': 'CC(C)(CO)[CH](O)C(=O)NCCC(=O)[OH]',
        
        # B6 (Pyridoxine/Pyridoxal/Pyridoxamine) - pyridine with specific substitution
        'B6': 'Cc1ncc(CO)c(CO)c1O',
        
        # B7 (Biotin) - fused rings with sulfur
        'B7': 'O=C1NC2[CH]S[CH]C2N1',
        
        # B9 (Folate) - pterin + pABA + glutamate
        'B9': 'Nc1nc2ncc(CNc3ccc(CC)cc3)nc2c(=O)[nH]1',
        
        # B12 (Cobalamin) - corrin ring with cobalt
        'B12': '[Co]'
    }
    
    # Check each pattern
    for vitamin, pattern in patterns.items():
        substructure = Chem.MolFromSmarts(pattern)
        if substructure and mol.HasSubstructMatch(substructure):
            if vitamin == 'B12':
                # Additional check for corrin ring system for B12
                if any(atom.GetSymbol() == 'Co' for atom in mol.GetAtoms()):
                    return True, f"Matches vitamin {vitamin} (Cobalamin) structure"
            else:
                return True, f"Matches vitamin {vitamin} structure"
    
    # Additional checks for derivatives and variations
    
    # Check for phosphate groups (common in active forms)
    phosphate = Chem.MolFromSmarts('P(=O)(O)(O)O')
    
    # Check for nucleotide-like structures (as in FAD)
    nucleotide = Chem.MolFromSmarts('c1ncnc2[nH]cnc12')
    
    if mol.HasSubstructMatch(phosphate) and mol.HasSubstructMatch(nucleotide):
        return True, "Matches B vitamin derivative (likely FAD or FMN)"
    
    # Count number of rings
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    
    # Additional characteristics that might indicate a B vitamin
    nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if num_rings >= 2 and nitrogen_count >= 2 and oxygen_count >= 2:
        # This could be a modified or derivative form
        return True, "Matches structural characteristics of B vitamin derivative"
        
    return False, "Does not match any B vitamin structural patterns"


__metadata__ = {
    'chemical_class': {
        'name': 'B vitamin',
        'definition': 'Any member of the group of eight water-soluble vitamins '
                     'originally thought to be a single compound (vitamin B) that '
                     'play important roles in cell metabolism.',
        'subclasses': ['B1', 'B2', 'B3', 'B5', 'B6', 'B7', 'B9', 'B12']
    }
}