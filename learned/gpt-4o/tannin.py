"""
Classifies: CHEBI:26848 tannin
"""
from rdkit import Chem

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    Tannins are characterized by complex polyphenolic structures with catechol/pyrogallol units 
    and sugar or sugar-like moieties, often forming glycosidic bonds or esterified forms.
    
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

    # General aromatic rings with hydroxyl groups
    aromatic_phenol = Chem.MolFromSmarts("a1ccc(O)cc1")
    aromatic_phenol_matches = mol.GetSubstructMatches(aromatic_phenol)
    if len(aromatic_phenol_matches) < 2:
        return False, "Less than two phenolic rings"

    # Check for catechol or pyrogallol moieties
    catechol_pattern = Chem.MolFromSmarts("c1c(O)cc(O)cc1")
    pyrogallol_pattern = Chem.MolFromSmarts("c1c(O)c(O)c(O)cc1")
    if not (mol.HasSubstructMatch(catechol_pattern) or mol.HasSubstructMatch(pyrogallol_pattern)):
        return False, "No catechol or pyrogallol found"

    # Check for the presence of sugar moieties or esterified forms
    potential_sugar_patterns = [
        Chem.MolFromSmarts("C1O[C@H](CO)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O"),  # Glucose-like and cyclic sugars
        Chem.MolFromSmarts("C1OC[C@@H](O)[C@@H](O)[C@H](O)[C@@H]1O"),  # Another sugar cyclic variant
        Chem.MolFromSmarts("C(=O)O[C@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O"),  # Esterified sugar
        Chem.MolFromSmarts("C(=O)[O|S|N]C")  # More general acyl or ester bonds
    ]
    
    sugar_found = any(mol.HasSubstructMatch(sugar) for sugar in potential_sugar_patterns)
    if not sugar_found:
        return False, "No identifiable sugar moieties found"

    return True, "Contains typical polyphenolic tannin structure with sugar moieties"