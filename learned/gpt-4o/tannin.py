"""
Classifies: CHEBI:26848 tannin
"""
from rdkit import Chem

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    A tannin is a complex polyphenolic compound often with catechol or pyrogallol residues and sugar moieties.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tannin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Search for multiple phenolic rings
    phenolic_pattern = Chem.MolFromSmarts("c1cc(O)ccc1")
    phenolic_matches = mol.GetSubstructMatches(phenolic_pattern)
    if len(phenolic_matches) < 2:
        return False, "Less than two phenolic rings"
    
    # Search for catechol or pyrogallol moieties
    catechol_pattern = Chem.MolFromSmarts("c1c(O)cc(O)cc1")
    pyrogallol_pattern = Chem.MolFromSmarts("c1c(O)c(O)c(O)cc1")
    if not (mol.HasSubstructMatch(catechol_pattern) or mol.HasSubstructMatch(pyrogallol_pattern)):
        return False, "No catechol or pyrogallol found"
    
    # Check for the presence of different sugar-like moieties
    sugar_patterns = [
        Chem.MolFromSmarts("C1OC(CO)C(O)C(O)C1O"),  # Glucose
        Chem.MolFromSmarts("C1C(O)C(O)C(O)C(O)C1O"),  # Other potential sugars
        Chem.MolFromSmarts("O[C@H]([C@@H]1O[C@@H]([C@H](O)[C@@H](O)[C@H]1O)C)C"),  # Glucoside-like connection
    ]
    
    sugar_found = any(mol.HasSubstructMatch(sugar) for sugar in sugar_patterns)
    if not sugar_found:
        return False, "No sugar connections found"
    
    return True, "Contains typical polyphenolic tannin structure with sugar moieties"