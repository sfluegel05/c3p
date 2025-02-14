"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
from rdkit import Chem

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.
    A nucleoside 5'-phosphate is a ribosyl or deoxyribosyl derivative of a pyrimidine
    or purine base with C-5 of the ribose ring being mono-, di-, tri- or tetra-phosphorylated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside 5'-phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a pattern for the sugar moiety (more generic version)
    sugar_pattern = Chem.MolFromSmarts("O[C@@H]1[C@H](CO)O[C@H](n)C1")

    # Pattern for 5'-phosphate group, could be one or more phosphate groups
    phosphate_pattern = Chem.MolFromSmarts("COP([O-])([O-])=O")
    
    # Check for the presence of the sugar moiety
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No recognized ribose or deoxyribose sugar detected"
    
    # Check for the presence of a phosphate group
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No recognizable phosphate group at the 5' position"

    # Check for nucleobase; purine or pyrimidine
    nucleobase_smarts = Chem.MolFromSmarts("n1cnc2c1ncnc2 |C|")

    if not mol.HasSubstructMatch(nucleobase_smarts):
        return False, "No recognizable nucleobase found"

    return True, "Contains ribose or deoxyribose sugar, nucleobase, and 5'-phosphate group"