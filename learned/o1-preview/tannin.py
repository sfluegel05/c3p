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
    
    # Calculate molecular weight - tannins typically have high molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight is {mol_wt:.1f} Da, which is less than 500 Da"
    
    # Count number of phenolic hydroxyl groups (aromatic hydroxyls)
    phenol_pattern = Chem.MolFromSmarts('c[OH]')
    phenol_matches = mol.GetSubstructMatches(phenol_pattern)
    num_phenols = len(phenol_matches)
    if num_phenols < 3:
        return False, f"Found {num_phenols} phenolic hydroxyl groups, need at least 3"
    
    # Count number of aromatic rings
    num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if num_aromatic_rings < 3:
        return False, f"Found {num_aromatic_rings} aromatic rings, need at least 3"
    
    # Check for ester linkages (common in tannins)
    ester_pattern = Chem.MolFromSmarts('C(=O)O')
    has_ester = mol.HasSubstructMatch(ester_pattern)
    
    # Check for glycosidic linkages (sugar moieties)
    glycosidic_pattern = Chem.MolFromSmarts('[OX2H][CX4H1]')
    has_glycosidic = mol.HasSubstructMatch(glycosidic_pattern)
    
    if not (has_ester or has_glycosidic):
        return False, "No ester or glycosidic linkages found"
    
    # All criteria met
    return True, "Contains multiple phenolic hydroxyls, aromatic rings, high molecular weight, and ester or glycosidic linkages"
    
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