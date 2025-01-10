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
    
    # Look for a broader nucleobase pattern (purines and pyrimidines with variations)
    purine_patt = Chem.MolFromSmarts('n1cnc2[nH]cnc12')  # Simplified for common purine bases
    pyrimidine_patt = Chem.MolFromSmarts('n1cncnc1')     # Simplified for pyrimidine bases

    nucleobase_found = mol.HasSubstructMatch(purine_patt) or mol.HasSubstructMatch(pyrimidine_patt)
    if not nucleobase_found:
        return False, "No nucleobase detected (common purine or pyrimidine)"

    # Look for ribose/deoxyribose-like structures, including flexibility for hydroxyl variations
    sugar_pattern = Chem.MolFromSmarts('C1O[C@@H](CO)[C@@H](O)[C@H]1O')  # Ribofuranose
    deoxy_sugar_pattern = Chem.MolFromSmarts('C1[C@H](O)[C@@H](O)[C@@H](O)C[O]1')  # Allow minor variations

    if not (mol.HasSubstructMatch(sugar_pattern) or mol.HasSubstructMatch(deoxy_sugar_pattern)):
        return False, "No ribose or deoxyribose sugar found or mismatched structure"

    # Check for phosphate groups including possible variations
    phosphate_patt = Chem.MolFromSmarts('P(=O)(O)O')  # Allow attachment to other oxygens or sugars
    if not mol.HasSubstructMatch(phosphate_patt):
        return False, "No phosphate group(s) detected"

    return True, "Molecule matches nucleoside phosphate structure due to presence of nucleobase, sugar, and phosphate"