"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate based on its SMILES string.
    A nucleoside phosphate has a nucleobase, a sugar moiety, and one or more phosphate groups attached.

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
    
    # Extended nucleobase patterns (covering purine and pyrimidine with common modifications)
    nucleobase_patterns = [
        Chem.MolFromSmarts("n1cnc2c1ncnc2"),  # Purine
        Chem.MolFromSmarts("n1c([nH])nc2c1[nH]cnc2"),  # Pyrimidine
        Chem.MolFromSmarts("c1ncnc2ncnc(N)c12")  # Adenine-like
        # Additional SMARTS for other nucleobase variants can be added
    ]
    
    nucleobase_found = any(mol.HasSubstructMatch(pattern) for pattern in nucleobase_patterns)
    if not nucleobase_found:
        return False, "No nucleobase found"
    
    # Refined sugar moiety detection
    sugar_patterns = [
        Chem.MolFromSmarts("C1OC(O)C(O)C(O)C1"),  # Ribose
        Chem.MolFromSmarts("C1OC(C)C(O)C1"),  # Deoxyribose, allowing some substitutions
        # Expand for other sugar forms as necessary
    ]
    
    sugar_found = any(mol.HasSubstructMatch(pattern) for pattern in sugar_patterns)
    if not sugar_found:
        return False, "No sugar moiety found"
    
    # Flexible phosphate group detection, allowing for mono- and polyphosphate
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)O")  # Single phosphate pattern
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"
    
    # Ensuring nucleobase and sugar and phosphate connectivity could involve
    # a graph-based connection check, but we will assume molecular subunits
    # provided in SMILES are connected in typical biological context

    return True, "Contains nucleobase, sugar moiety, and phosphate group"