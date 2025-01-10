"""
Classifies: CHEBI:50699 oligosaccharide
"""
"""
Classifies: oligosaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdqueries

def is_oligosaccharide(smiles: str):
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string.
    An oligosaccharide is composed of monosaccharide units joined by glycosidic linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oligosaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Kekulize the molecule (helps with certain valence issues)
    try:
        Chem.Kekulize(mol)
    except Chem.KekulizeException as e:
        pass  # Kekulization may fail for some molecules, but we can continue

    # Define monosaccharide pattern (5 or 6 membered ring with oxygen and hydroxyl groups)
    monosaccharide_patterns = [
        Chem.MolFromSmarts('[C;H1,H2][O][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2]'),  # Pyranose ring
        Chem.MolFromSmarts('[C;H1,H2][O][C;H1,H2][C;H1,H2][C;H1,H2]')             # Furanose ring
    ]
    
    # Find monosaccharide units
    monosaccharide_matches = []
    for pattern in monosaccharide_patterns:
        matches = mol.GetSubstructMatches(pattern)
        monosaccharide_matches.extend(matches)
    
    num_monosaccharides = len(set([tuple(sorted(match)) for match in monosaccharide_matches]))
    
    if num_monosaccharides < 2:
        return False, f"Found {num_monosaccharides} monosaccharide unit(s), need at least 2 for oligosaccharide"
    
    # Identify glycosidic linkages (C-O-C between rings)
    glycosidic_bond_pattern = Chem.MolFromSmarts('[C;R][O;$(*-[#6]);R][$([C;R]),$([C;r3,r4,r5,r6])]')
    glycosidic_bonds = mol.GetSubstructMatches(glycosidic_bond_pattern)
    
    num_glycosidic_bonds = len(glycosidic_bonds)
    
    if num_glycosidic_bonds < 1:
        return False, "No glycosidic linkages found between monosaccharide units"
    
    # Check if the number of glycosidic bonds corresponds to number of monosaccharides minus one
    if num_glycosidic_bonds < num_monosaccharides - 1:
        return False, f"Insufficient glycosidic linkages: found {num_glycosidic_bonds}, expected at least {num_monosaccharides - 1}"
    
    # Additional check: ensure that the molecule is not too large (exclude polysaccharides)
    if num_monosaccharides > 10:
        return False, f"Found {num_monosaccharides} monosaccharide units, which may indicate a polysaccharide"
    
    return True, f"Contains {num_monosaccharides} monosaccharide units connected via glycosidic linkages"

__metadata__ = {   'chemical_class': {   'id': None,
                              'name': 'oligosaccharide',
                              'definition': 'A compound in which monosaccharide units are joined by glycosidic linkages. The term is commonly used to refer to a defined structure as opposed to a polymer of unspecified length or a homologous mixture. When the linkages are of other types the compounds are regarded as oligosaccharide analogues.',
                              'parents': []},
        'config': {   'llm_model_name': 'lbl/claude-sonnet',
                      'f1_threshold': 0.8,
                      'max_attempts': 5,
                      'max_positive_instances': None,
                      'max_positive_to_test': None,
                      'max_negative_to_test': None,
                      'max_positive_in_prompt': 50,
                      'max_negative_in_prompt': 20,
                      'max_instances_in_prompt': 100,
                      'test_proportion': 0.1},
        'message': None,
        'attempt': 0,
        'success': None,
        'best': None,
        'error': '',
        'stdout': None,
        'num_true_positives': None,
        'num_false_positives': None,
        'num_true_negatives': None,
        'num_false_negatives': None,
        'num_negatives': None,
        'precision': None,
        'recall': None,
        'f1': None,
        'accuracy': None}