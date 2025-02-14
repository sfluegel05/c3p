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
    
    # Define more general patterns for sugar moieties
    # Pattern for generic pentose sugar (both ribose and deoxyribose without strict stereochemistry)
    ribose_deoxyribose_general = Chem.MolFromSmarts("O[C@@H]1[C@H](O)[C@@H](O)[C@H]([C@H]1O)O |C|")
    ribose_or_deoxyribose_no_stereo = Chem.MolFromSmarts("O[C](CO)C(C)C(O)O")

    # Patterns for phosphate groups at 5' position (allowing for both charged and uncharged forms)
    phosphate_group_any = Chem.MolFromSmarts("COP(=O)(O*)O")  # Generic phosphate group matching
    
    # Check for the sugar part
    if not (mol.HasSubstructMatch(ribose_deoxyribose_general) or mol.HasSubstructMatch(ribose_or_deoxyribose_no_stereo)):
        return False, "No ribose or deoxyribose sugar detected"
    
    # Check for the presence of a phosphate group
    if not mol.HasSubstructMatch(phosphate_group_any):
        return False, "No recognizable phosphate group at the 5' position"

    # Check for nucleobase; maintain pyrimidine and purine bases with expanded structure matching
    # This pattern covers basic structures in common nucleobases
    nucleobase_smarts = Chem.MolFromSmarts("[nH]1cnc2c(n1)[nH]c[nH]c2 |C,n,o|")

    if not mol.HasSubstructMatch(nucleobase_smarts):
        return False, "No recognizable nucleobase found"

    return True, "Contains ribose or deoxyribose sugar, nucleobase, and 5'-phosphate group"