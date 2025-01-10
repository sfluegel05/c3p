"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
"""
Classifies: pyrimidine deoxyribonucleoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_pyrimidine_deoxyribonucleoside(smiles: str):
    """
    Determines if a molecule is a pyrimidine deoxyribonucleoside based on its SMILES string.
    A pyrimidine deoxyribonucleoside consists of a pyrimidine base connected to a 
    deoxyribose sugar via an N-glycosidic bond.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a pyrimidine deoxyribonucleoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for deoxyribose sugar pattern
    # More flexible pattern that matches 5-membered sugar ring
    # with a 2' deoxy position
    sugar_pattern = Chem.MolFromSmarts("[O;R1]1[CH2][CH1][CH1][CH1]1")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No deoxyribose sugar ring found"

    # Check for pyrimidine base patterns
    # Match common pyrimidine bases (uracil, thymine, cytosine)
    uracil_pattern = Chem.MolFromSmarts("O=c1[nH]c(=O)ccn1")
    thymine_pattern = Chem.MolFromSmarts("O=c1[nH]c(=O)c(C)cn1")
    cytosine_pattern = Chem.MolFromSmarts("[NH2]c1nc(=O)ccn1")
    
    has_pyrimidine = any(mol.HasSubstructMatch(pattern) for pattern in 
                        [uracil_pattern, thymine_pattern, cytosine_pattern])
    
    if not has_pyrimidine:
        # Fallback to more general pyrimidine pattern
        general_pyrimidine = Chem.MolFromSmarts("n1c[c,n]c[c,n]c1")
        if not mol.HasSubstructMatch(general_pyrimidine):
            return False, "No pyrimidine base found"

    # Check for N-glycosidic bond between sugar and base
    # More flexible pattern that allows for various substitutions
    glycosidic_pattern = Chem.MolFromSmarts("[O;R1]1[CH2][CH1][CH1][CH1]1[N;R1]2[#6][#6][#6][#6][#6]2")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No N-glycosidic bond between sugar and base"

    # Check for primary alcohol (5' position) - allow for modifications
    primary_oh_pattern = Chem.MolFromSmarts("[CH2][O,N,S]")
    if not mol.HasSubstructMatch(primary_oh_pattern):
        return False, "Missing 5' substituent"

    # Check for 3' position - allow for various substituents
    three_prime_pattern = Chem.MolFromSmarts("[CH1;R1][OH1,O,N,F]")
    if not mol.HasSubstructMatch(three_prime_pattern):
        return False, "Missing 3' substituent"

    return True, "Contains pyrimidine base connected to deoxyribose via N-glycosidic bond"