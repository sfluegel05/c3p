"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
from rdkit import Chem

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.
    A nucleoside 5'-phosphate consists of a ribosyl or deoxyribosyl sugar linked to
    a nucleobase with a phosphate group at the 5' position of the sugar.

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
    
    # Define SMARTS patterns for sugar moieties
    # Pattern for generic pentose sugar (catch both ribose and deoxyribose with stereochemistry)
    ribose_deoxyribose_pattern = Chem.MolFromSmarts("[C@H]1(O)[C@@H](O)[C@H](O)C[C@@H]1O")

    # Patterns for phosphate group at 5' position
    phosphate_pattern_mono = Chem.MolFromSmarts("COP([O-])(=O)[O-]")
    phosphate_pattern_di = Chem.MolFromSmarts("COP([O-])(=O)OP([O-])(=O)[O-]")
    
    # Check for the sugar part
    if not mol.HasSubstructMatch(ribose_deoxyribose_pattern):
        return False, "No ribose or deoxyribose sugar detected"
    
    # Check for the presence of a phosphate group
    if not (mol.HasSubstructMatch(phosphate_pattern_mono) or mol.HasSubstructMatch(phosphate_pattern_di)):
        return False, "No phosphate group at the 5' position"

    # Simplified nucleobase presence check using generic pyrimidine and purine base structures
    pyrimidine_base = Chem.MolFromSmarts("c1cncnc1")
    purine_base = Chem.MolFromSmarts("n1cnc2c1ncnc2")
    
    if not (mol.HasSubstructMatch(pyrimidine_base) or mol.HasSubstructMatch(purine_base)):
        return False, "No recognizable nucleobase found"

    return True, "Contains ribose or deoxyribose sugar, nucleobase, and 5'-phosphate group"