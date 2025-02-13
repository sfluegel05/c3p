"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
from rdkit import Chem

def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate based on its SMILES string.
    A nucleoside phosphate is characterized by a nucleobase attached to a sugar with 
    one or more phosphate groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a nucleoside phosphate, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Comprehensive list of nucleobase patterns (SMARTS)
    nucleobases = [
        'c1nc[nH]c2c1ncnc2N',       # Adenine
        'Nc1ncnc2n(ccn(c12))',      # Guanine 
        'c1cc(=O)[nH]c(n1)[C@H]2O', # Thymine
        'O=C1NC=NC2=C1[N]C[C@H]2O', # Uracil
        'n1cc(c(=O)[nH]c1)[C@@H]2O' # Cytosine
    ]

    # Check for nucleobase structure in the molecule
    nucleobase_detected = False
    for base in nucleobases:
        base_pattern = Chem.MolFromSmarts(base)
        if mol.HasSubstructMatch(base_pattern):
            nucleobase_detected = True
            break

    if not nucleobase_detected:
        return False, "Nucleobase not found in molecule"

    # Phosphate group patterns; including polyphosphates and unique linkages
    phosphate_patterns = [
        '[OP](=O)(O)[O-]',      # Mono phosphate ion
        '[O-]P(=O)(O)O[C@@H]',  # Additional linkage allowing for phosphate esters
    ]
    
    # Check if at least one phosphate group is present
    phosphate_detected = False
    for phosphate in phosphate_patterns:
        phosphate_pattern = Chem.MolFromSmarts(phosphate)
        if mol.HasSubstructMatch(phosphate_pattern):
            phosphate_detected = True
            break

    if not phosphate_detected:
        return False, "No phosphate group found in the molecule"

    # Check for sugar connection (e.g., ribose or deoxyribose)
    sugar_patterns = [
        '[C@H]1(O)C(O)C(O)C(O1)',   # Ribose
        '[C@H]1(O)C[C@H](O)C(O1)O'  # Deoxyribose
    ]
    
    sugar_detected = False
    for sugar in sugar_patterns:
        sugar_pattern = Chem.MolFromSmarts(sugar)
        if mol.HasSubstructMatch(sugar_pattern):
            sugar_detected = True
            break

    if not sugar_detected:
        return False, "No sugar structure found in nucleoside phosphate"

    return True, "Contains a nucleobase, sugar, and phosphate group, classified as nucleoside phosphate"