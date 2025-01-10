"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
from rdkit import Chem

def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate based on its SMILES string.
    
    A nucleoside phosphate is characterized by a nucleobase connected to a sugar
    moiety, which is phosphorylated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define more flexible patterns
    # Broad nucleobase pattern for common bases and their derivatives
    purine_patterns = [
        Chem.MolFromSmarts('n1cnc2ncnc12'),  # Common purine
        Chem.MolFromSmarts('c1ncnc2[nH]cnc12'),  # Modified purine
    ]
    pyrimidine_patterns = [
        Chem.MolFromSmarts('n1c([nH])ccnc1'),  # Common pyrimidine
        Chem.MolFromSmarts('n1cncn[nH]1'),    # Modified pyrimidine
    ]

    nucleobase_found = any(mol.HasSubstructMatch(pat) for pat in purine_patterns + pyrimidine_patterns)
    if not nucleobase_found:
        return False, "No nucleobase detected (general purine or pyrimidine)"

    # Check for sugar moiety, more flexible and allowing for common modifications
    sugar_patterns = [
        Chem.MolFromSmarts('C1OCC(O)C1O'),  # Ribose/deoxyribose core without strict stereochemistry
        Chem.MolFromSmarts('OCC1OC(CO)C(O)C1'),  # Other sugar forms
    ]

    sugar_found = any(mol.HasSubstructMatch(pat) for pat in sugar_patterns)
    if not sugar_found:
        return False, "No appropriate sugar detected (ribose variants)"

    # Look for phosphate groups (allow different attachments)
    phosphate_patterns = [
        Chem.MolFromSmarts('P(=O)(O)(O)'),  # Mono-phosphate
        Chem.MolFromSmarts('OP(=O)(O)OP(=O)(O)O'),  # Di-phosphate
        Chem.MolFromSmarts('OP(=O)(O)OP(=O)(O)OP(=O)(O)O'),  # Tri-phosphate
    ]

    phosphate_found = any(mol.HasSubstructMatch(pat) for pat in phosphate_patterns)
    if not phosphate_found:
        return False, "No phosphate group(s) detected"

    return True, "Molecule matches nucleoside phosphate structure due to presence of nucleobase, sugar, and phosphate"