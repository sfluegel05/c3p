"""
Classifies: CHEBI:26848 tannin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    Tannins are polyphenolic compounds with multiple hydroxyl groups on aromatic rings,
    often with glycosidic or ester linkages.

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

    # Count aromatic rings
    aromatic_rings = sum(1 for ring in Chem.GetSymmSSSR(mol) if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring))
    if aromatic_rings < 1:
        return False, "No aromatic rings found"

    # Count phenolic hydroxyl groups (-OH attached to aromatic rings)
    phenolic_oh_pattern = Chem.MolFromSmarts("[OH][c]")
    phenolic_oh_matches = mol.GetSubstructMatches(phenolic_oh_pattern)
    if len(phenolic_oh_matches) < 2:
        return False, f"Found {len(phenolic_oh_matches)} phenolic hydroxyl groups, need at least 2"

    # Check for glycosidic or ester linkages
    glycosidic_pattern = Chem.MolFromSmarts("[OX2][C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O1)")  # Glycosidic linkage
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")  # Ester linkage
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(glycosidic_matches) == 0 and len(ester_matches) == 0:
        return False, "No glycosidic or ester linkages found"

    # Check molecular weight - tannins are typically large molecules
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for tannin"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 10:
        return False, "Too few carbons for tannin"
    if o_count < 5:
        return False, "Too few oxygens for tannin"

    return True, "Contains multiple phenolic hydroxyl groups on aromatic rings with glycosidic or ester linkages"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33248',
                          'name': 'tannin',
                          'definition': 'Any of a group of astringent polyphenolic vegetable principles or compounds, chiefly complex glucosides of catechol and pyrogallol.',
                          'parents': ['CHEBI:33248', 'CHEBI:33248']},
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
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199}