"""
Classifies: CHEBI:18241 2'-deoxyribonucleoside 5'-monophosphate
"""
"""
Classifies: CHEBI:73333 2'-deoxyribonucleoside 5'-monophosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2__deoxyribonucleoside_5__monophosphate(smiles: str):
    """
    Determines if a molecule is a 2'-deoxyribonucleoside 5'-monophosphate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 2'-deoxyribonucleoside 5'-monophosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define key structural patterns
    patterns = {
        # Deoxyribose sugar with 5' phosphate connection
        # Matches furanose ring with 3'-OH, no 2'-OH, and 5'-O-P
        'deoxyribose_5p': '[O,OH1]-[P](=[O])([O,OH1])[O,OH1]-[CH2]-[CH]1-[O]-[CH](-[C,N])-[CH](-[OH1])-[CH2]1',
        
        # Alternative pattern for charged phosphate forms
        'deoxyribose_5p_charged': '[O-]-[P](=[O])([O,O-])[O]-[CH2]-[CH]1-[O]-[CH](-[C,N])-[CH](-[OH1])-[CH2]1',
        
        # Basic nucleobase patterns - more general
        'base_pattern': '[C,N]1[C,N][C,N][C,N][C,N]1',  # General ring pattern for nucleobases
    }
    
    # Convert patterns to RDKit molecules
    smart_patterns = {}
    for name, pattern in patterns.items():
        smart_mol = Chem.MolFromSmarts(pattern)
        if smart_mol is None:
            return False, f"Invalid SMARTS pattern for {name}"
        smart_patterns[name] = smart_mol
    
    # Check for deoxyribose with 5' phosphate
    has_deoxy_5p = mol.HasSubstructMatch(smart_patterns['deoxyribose_5p']) or \
                   mol.HasSubstructMatch(smart_patterns['deoxyribose_5p_charged'])
    
    if not has_deoxy_5p:
        return False, "No deoxyribose sugar with 5' phosphate found"
    
    # Check for nucleobase
    if not mol.HasSubstructMatch(smart_patterns['base_pattern']):
        return False, "No nucleobase ring system found"
    
    # Count phosphorus atoms - should only have one
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count != 1:
        return False, f"Found {p_count} phosphorus atoms, should be exactly 1"
    
    # Additional validation checks
    
    # Check for reasonable number of oxygen atoms (at least 5: 3 from phosphate, 2 from sugar)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 5:
        return False, "Insufficient oxygen atoms for structure"
    
    # Check for reasonable molecular size
    if mol.GetNumAtoms() < 15:
        return False, "Molecule too small to be a nucleotide"
    
    return True, "Contains deoxyribose sugar with 5' phosphate and nucleobase"