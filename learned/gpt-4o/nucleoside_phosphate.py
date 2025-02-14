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

    # Phosphate group patterns
    phosphate_patterns = [
        '[OX1P](=O)(O)[O]',      # Ester-linked phosphate (covers both mono- and poly-phosphate)
        '[OX1P](=O)(O)(O)',      # Phosphate group attached, including a potential -OP linkage
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

    # Check for sugar connection (e.g., ribose or deoxyribose typical in nucleotides)
    sugar_pattern = Chem.MolFromSmarts('C1C(O)C(O)C(O)C1O')  # Simple pattern for ribose/deoxyribose
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar structure found in nucleoside phosphate"

    return True, "Contains a nucleobase, sugar, and phosphate group, classified as nucleoside phosphate"