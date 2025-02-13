"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate based on its SMILES string.
    A nucleoside phosphate is defined as a nucleobase-containing molecular entity that is 
    a nucleoside in which one or more of the sugar hydroxy groups has been converted into 
    a mono- or poly-phosphate.

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

    # List of nucleobase patterns (SMARTS)
    nucleobases = [
        '[nH]1cnc2c1ncnc2N',          # Adenine
        'n1ccn(c(=O)[nH]1)[C@H]2O',   # Cytosine
        'Nc1nc2n(cnc2c(=O)[nH]1)',    # Guanine 
        'c1cc(=O)nc(n1)[C@H]2O',      # Thymine
        'O=C1NC=NC2=C1C@H(O)C@H>2)O'  # Uracil
    ]

    # Check for nucleobase structure
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
        'P(=O)(O)(O)',        # monophosphate
        'OP(=O)(O)[O-]',      # deprotonated monophosphate variant
        'P(=O)(O)(O)O',       # monoester phosphate
        'P(=O)(O)(O)OCC'      # phosphate linked to sugar
    ]
    
    # Check if at least one phosphate group is present
    phosphate_detected = False
    for phosphate in phosphate_patterns:
        phosphate_pattern = Chem.MolFromSmarts(phosphate)
        if mol.HasSubstructMatch(phosphate_pattern):
            phosphate_detected = True
            break

    if not phosphate_detected:
        return False, "No phosphate group found in molecule"

    return True, "Contains a nucleobase and phosphate group, thus classified as nucleoside phosphate"