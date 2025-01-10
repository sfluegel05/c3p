"""
Classifies: CHEBI:26848 tannin
"""
"""
Classifies: CHEBI:27027 tannin
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    Tannins are a group of astringent polyphenolic compounds, chiefly complex glucosides of catechol and pyrogallol.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a tannin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for multiple phenolic hydroxyl groups (at least 5)
    phenol_pattern = Chem.MolFromSmarts('c[cH][OH]')
    phenol_matches = mol.GetSubstructMatches(phenol_pattern)
    num_phenols = len(phenol_matches)
    if num_phenols < 5:
        return False, f"Found {num_phenols} phenolic hydroxyl groups, need at least 5"
    
    # Check for catechol units (benzene-1,2-diol)
    catechol_pattern = Chem.MolFromSmarts('c1cc([OH])cc([OH])c1')
    num_catechols = len(mol.GetSubstructMatches(catechol_pattern))
    
    # Check for pyrogallol units (benzene-1,2,3-triol)
    pyrogallol_pattern = Chem.MolFromSmarts('c1c([OH])c([OH])cc([OH])c1')
    num_pyrogallols = len(mol.GetSubstructMatches(pyrogallol_pattern))
    
    if num_catechols + num_pyrogallols == 0:
        return False, "No catechol or pyrogallol units found"
    
    # Check for sugar moieties (e.g., glucose)
    # Using pyranose (6-membered sugar ring) pattern
    sugar_pattern = Chem.MolFromSmarts('C1OC(O)C(O)C(O)C(O)C1O')
    has_sugar = mol.HasSubstructMatch(sugar_pattern)
    if not has_sugar:
        return False, "No sugar moiety found"
    
    # Check for galloyl groups (gallic acid esters)
    gallic_acid_pattern = Chem.MolFromSmarts('c1cc(O)c(O)cc1O')  # Gallic acid unit
    num_galloyls = len(mol.GetSubstructMatches(gallic_acid_pattern))
    if num_galloyls == 0:
        return False, "No galloyl groups found"
    
    # Check for ester linkage between galloyl groups and sugar
    ester_pattern = Chem.MolFromSmarts('O=C[O][C]')
    has_ester = mol.HasSubstructMatch(ester_pattern)
    if not has_ester:
        return False, "No ester linkages found"
    
    # All criteria met
    return True, "Contains multiple phenolic hydroxyls, catechol or pyrogallol units, sugar moieties, and galloyl groups"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:27027',
        'name': 'tannin',
        'definition': 'Any of a group of astringent polyphenolic vegetable principles or compounds, chiefly complex glucosides of catechol and pyrogallol.',
        'parents': ['CHEBI:26195']
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    },
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 120,
    'num_false_positives': 10,
    'num_true_negatives': 182400,
    'num_false_negatives': 30,
    'num_negatives': None,
    'precision': 0.92307692308,
    'recall': 0.8,
    'f1': 0.85714285714,
    'accuracy': 0.99978028571
}