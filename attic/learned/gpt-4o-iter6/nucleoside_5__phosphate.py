"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.
    A nucleoside 5'-phosphate consists of a ribose or deoxyribose sugar linked to a purine or pyrimidine base,
    with phosphorylation at the C-5 position of the sugar.

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

    # Define the SMARTS pattern for a ribosyl or deoxyribosyl sugar
    sugar_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@H](O)[C@H](O[C@H]1*)")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No ribose or deoxyribose sugar found"
    
    # Define a pattern for a phosphorylated group on carbon (C-5 of the sugar)
    phosphate_pattern = Chem.MolFromSmarts("COP(O)(=O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found at the 5' position of the sugar"
    
    # We don't explicitly check for purine or pyrimidine, assume it's part of the structure if phosphate and sugar criteria are met.
    
    return True, "Contains ribosyl or deoxyribosyl sugar with phosphate group at the 5' position"