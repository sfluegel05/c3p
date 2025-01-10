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
    
    # Define SMARTS patterns
    # Catechol unit: benzene ring with adjacent hydroxyls at positions 1 and 2
    catechol_pattern = Chem.MolFromSmarts('c1c(O)cccc1O')
    # Pyrogallol unit: benzene ring with hydroxyls at positions 1,2,3
    pyrogallol_pattern = Chem.MolFromSmarts('c1c(O)c(O)c(O)cc1')
    # Galloyl ester: ester of gallic acid
    galloyl_ester_pattern = Chem.MolFromSmarts('O=C(O)c1c(O)ccc(O)c1')
    # Glycosidic linkage: sugar linkage (C-O-C between carbons)
    glycosidic_pattern = Chem.MolFromSmarts('[CX4][OX2H][CX4]')
    # Flavonoid unit: approximate pattern for flavan-3-ol core
    flavonoid_pattern = Chem.MolFromSmarts('C1[C@@H](O)[C@@H](C2=CC=CC=C2)OC3=CC=CC=C13')
    
    # Count catechol units
    num_catechol = len(mol.GetSubstructMatches(catechol_pattern))
    # Count pyrogallol units
    num_pyrogallol = len(mol.GetSubstructMatches(pyrogallol_pattern))
    # Count galloyl ester groups
    num_galloyl = len(mol.GetSubstructMatches(galloyl_ester_pattern))
    # Count glycosidic linkages
    num_glycosidic = len(mol.GetSubstructMatches(glycosidic_pattern))
    # Count flavonoid units
    num_flavonoid = len(mol.GetSubstructMatches(flavonoid_pattern))
    
    # Total number of phenolic units
    num_phenolic_units = num_catechol + num_pyrogallol + num_galloyl + num_flavonoid
    
    # Check if molecule has multiple phenolic units
    if num_phenolic_units < 2:
        return False, f"Found {num_phenolic_units} phenolic units, need at least 2"
    
    # Tannins usually have multiple galloyl esters or flavonoid units
    if num_galloyl + num_flavonoid < 2:
        return False, f"Found {num_galloyl} galloyl esters and {num_flavonoid} flavonoid units, need at least 2 combined"
    
    # Check for glycosidic linkage
    if num_glycosidic < 1:
        return False, "No glycosidic linkages found"
    
    # Verify high degree of polymerization
    num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if num_aromatic_rings < 4:
        return False, f"Found {num_aromatic_rings} aromatic rings, need at least 4 for tannins"
    
    # All criteria met
    return True, "Contains multiple phenolic units, galloyl esters/flavonoid units, glycosidic linkages, and high degree of polymerization"

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
    }
}