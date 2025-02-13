"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.
    Requires presence of a purine or pyrimidine base and a ribose or deoxyribose
    attached to a phosphate group at the 5' position.

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

    # Identify sugar component (ribose or deoxyribose)
    ribose_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H](C)O[C@H]([C@@H]1O)O")
    deoxyribose_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H](C)O[C@H](C1)O")
    if not mol.HasSubstructMatch(ribose_pattern) and not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No ribose or deoxyribose sugar found"

    # Identify purine or pyrimidine base
    purine_pattern = Chem.MolFromSmarts("c1ncnc2n(cnc12)")
    pyrimidine_pattern = Chem.MolFromSmarts("c1ncnc2n(cnc12)")
    if not (mol.HasSubstructMatch(purine_pattern) or mol.HasSubstructMatch(pyrimidine_pattern)):
        return False, "No purine or pyrimidine base found"

    # Identify phosphate group at C-5'
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 1:
        return False, "No phosphate group found at C-5'"

    return True, "Contains a ribosyl or deoxyribosyl derivative of a purine/pyrimidine base with a 5'-phosphate"