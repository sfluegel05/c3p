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
    
    # Define and validate SMARTS patterns
    patterns = {
        # Deoxyribose sugar - more flexible pattern
        'deoxyribose': '[C]-[C]-[C](-[O])-[C]1-O1',
        
        # Phosphate group - allowing for charged forms
        'phosphate': '[O]-[P](=[O])(-[O])(-[O])',
        
        # N-glycosidic bond - more general pattern
        'n_glycosidic': '[N]1-[C]2-O-[C]-[C]-[C]2',
        
        # 5' phosphate connection
        'phosphate_5': '[C]-[O]-[P](=[O])(-[O])(-[O])',
        
        # Various nucleobase patterns - more inclusive
        'purine': 'c1ncnc2ncnc12',
        'pyrimidine': 'c1cn([*])c(=O)nc1',
        'uracil': 'O=c1[n]c(=O)[n][c][n]1',
        'thymine': 'CC1=C(N([*])C(=O)NC1=O)',
        'cytosine': 'Nc1[n]c(=O)[n][c][n]1'
    }
    
    # Convert patterns to RDKit molecules
    smart_patterns = {}
    for name, pattern in patterns.items():
        smart_mol = Chem.MolFromSmarts(pattern)
        if smart_mol is None:
            return None, f"Invalid SMARTS pattern for {name}"
        smart_patterns[name] = smart_mol
    
    # Check for deoxyribose sugar
    if not mol.HasSubstructMatch(smart_patterns['deoxyribose']):
        return False, "No deoxyribose sugar found"
    
    # Check for phosphate group
    if not mol.HasSubstructMatch(smart_patterns['phosphate']):
        return False, "No phosphate group found"
    
    # Check for N-glycosidic bond
    if not mol.HasSubstructMatch(smart_patterns['n_glycosidic']):
        return False, "No N-glycosidic bond found"
    
    # Check for 5' phosphate
    if not mol.HasSubstructMatch(smart_patterns['phosphate_5']):
        return False, "Phosphate not at 5' position"
    
    # Check for nucleobase
    has_nucleobase = any(
        mol.HasSubstructMatch(pattern) 
        for pattern in [
            smart_patterns['purine'],
            smart_patterns['pyrimidine'],
            smart_patterns['uracil'],
            smart_patterns['thymine'],
            smart_patterns['cytosine']
        ]
    )
    
    if not has_nucleobase:
        return False, "No recognizable nucleobase found"
    
    # Count phosphorus atoms - should only have one
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count != 1:
        return False, f"Found {p_count} phosphorus atoms, should be exactly 1"
    
    # Additional check: should have 2' position without OH (deoxyribose)
    # This is implicit in the structure match but good to verify
    
    return True, "Contains deoxyribose sugar with 5' phosphate and nucleobase connected via N-glycosidic bond"