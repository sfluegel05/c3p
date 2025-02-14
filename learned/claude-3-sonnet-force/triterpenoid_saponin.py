"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
"""
Classifies: CHEBI:36624 triterpenoid saponin
A terpene glycoside in which the terpene moiety is a triterpenoid.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid saponin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for triterpenoid backbone patterns
    triterpene_backbones = [
        Chem.MolFromSmarts("[C@@]12[C@H]([C@@H]3[C@]([C@@H]4[C@H](C[C@@H]5[C@@]4(C)C[C@@H](O)[C@@]56C)C3(C)C)C2)C[C@@H]([C@@]1(C)C)O",  # oleanane
        Chem.MolFromSmarts("[C@@]12[C@H]([C@@H]3[C@]([C@@H]4[C@H](C[C@@H]5[C@@]4(C)C[C@@H](O)[C@@]56C)C3(C)C)C2)C[C@@H]([C@@]1(C)C)O",  # ursane
        Chem.MolFromSmarts("[C@@]12[C@H]([C@@H]3[C@]([C@@H]4[C@H](C[C@@H]5[C@@]4(C)C[C@@H](O)[C@@]56C)C3(C)C)C2)C[C@@H]([C@@]1(C)C)O"   # lupane
    ]
    has_triterpene_backbone = any(mol.HasSubstructMatch(pattern) for pattern in triterpene_backbones)
    if not has_triterpene_backbone:
        return False, "No triterpenoid backbone found"

    # Look for sugar moieties
    sugar_pattern = Chem.MolFromSmarts("[OX2r3]")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    has_sugar_moiety = len(sugar_matches) > 0
    if not has_sugar_moiety:
        return False, "No sugar moiety found"

    # Check for long chains (potential sugar chains)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 6:
        return False, "Potential sugar chains too short"

    # Check molecular weight - triterpenoid saponins typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500 or mol_wt > 1500:
        return False, "Molecular weight outside typical range for triterpenoid saponins"

    # Check elemental composition - should contain C, H, O
    atoms = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    has_c = 6 in atoms
    has_h = 1 in atoms
    has_o = 8 in atoms
    if not (has_c and has_h and has_o):
        return False, "Missing expected elements (C, H, O)"

    return True, "Contains triterpenoid backbone and sugar moiety"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:36624',
        'name': 'triterpenoid saponin',
        'definition': 'A terpene glycoside in which the terpene moiety is a triterpenoid.',
        'parents': ['CHEBI:35701', 'CHEBI:36623']
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
    'attempt': 5,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 149,
    'num_false_positives': 0,
    'num_true_negatives': 182432,
    'num_false_negatives': 3,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 0.9801980198019802,
    'f1': 0.9900498508175089,
    'accuracy': 0.9999817424242424
}