"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate based on its SMILES string.
    Nucleoside phosphates contain a nitrogenous base bonded to a sugar, with one or more phosphate groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a nucleoside phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string into RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Updated SMARTS pattern to detect phosphate groups (including di- and triphosphates)
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    polyphosphate_pattern = Chem.MolFromSmarts("P(=O)(O)OP(=O)(O)O")

    if not (mol.HasSubstructMatch(phosphate_pattern) or mol.HasSubstructMatch(polyphosphate_pattern)):
        return False, "No phosphate group or polyphosphate group found"

    # Broadening the pattern to identify pentose sugars with possible variations
    sugar_patterns = [
        Chem.MolFromSmarts("C1[C@@H](O)[C@H](O)[C@H](O)[C@@H]1O"),  # Ribose or Deoxyribose
        Chem.MolFromSmarts("[C@H]1(O)[C@@H](O)[C@H](O)[C@@H]([C@H]1O)"),  # Alt form
    ]

    if not any(mol.HasSubstructMatch(sp) for sp in sugar_patterns):
        return False, "No recognizable sugar moiety (ribose or similar) found"
    
    # Expanded base detection, including potential modified bases
    base_patterns = [
        Chem.MolFromSmarts("c1ncnc2ncnc12"),  # Purine structure
        Chem.MolFromSmarts("c1[nH]cnc2[nH]c[nH]c12"),  # Another form of specific pyrimidine-like pattern
        Chem.MolFromSmarts("c1cncnc1"),  # Simplified pyridine/pyrimidine
    ]
    
    if not any(mol.HasSubstructMatch(bp) for bp in base_patterns):
        return False, "No recognizable nitrogenous base found"
    
    return True, "Molecule contains features consistent with a nucleoside phosphate"

# Example test case
test_smiles = "Nc1ccn([C@@H]2O[C@H](CO)[C@@H](OP(O)(O)=O)[C@H]2O)c(=O)n1"
result, reason = is_nucleoside_phosphate(test_smiles)
print(result, reason)