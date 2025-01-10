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
    
    # Look for nucleobase patterns (purines and pyrimidines)
    purine_pattern = Chem.MolFromSmarts("n1cnc2c1ncnc2")  # Common purine structure
    pyrimidine_pattern = Chem.MolFromSmarts("n1c([nH])nc2c1[nH]cnc2")  # Common pyrimidine structure
    if not (mol.HasSubstructMatch(purine_pattern) or mol.HasSubstructMatch(pyrimidine_pattern)):
        return False, "No nucleobase found"
    
    # Look for sugar moiety patterns (ribose or deoxyribose)
    sugar_pattern = Chem.MolFromSmarts("C1OC(O)C(O)C1")  # Simple sugar ring
    deoxysugar_pattern = Chem.MolFromSmarts("C1OC1")  # Simple deoxysugar ring
    if not (mol.HasSubstructMatch(sugar_pattern) or mol.HasSubstructMatch(deoxysugar_pattern)):
        return False, "No sugar moiety found"
    
    # Look for phosphate groups (PO4)
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")  # Phosphate group
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) == 0:
        return False, "No phosphate group found"
    
    # Check for connection between nucleobase, sugar, and phosphate
    # For simplicity, assume if each part is present, they are connected correctly
    # Advanced checks would require analyzing connectivity further

    return True, "Contains nucleobase, sugar moiety, and phosphate group"