"""
Classifies: CHEBI:18379 nitrile
"""
"""
Classifies: CHEBI:17116 nitrile
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_nitrile(smiles: str):
    """
    Determines if a molecule is a nitrile based on its SMILES string.
    A nitrile is a compound having the structure RC#N, where the nitrile group (-C≡N) is attached to a carbon atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrile, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for nitrile pattern (C-C#N)
    nitrile_pattern = Chem.MolFromSmarts("[C;!$(C=*)]#N")
    if not mol.HasSubstructMatch(nitrile_pattern):
        return False, "No nitrile group (C-C#N) found"
    
    # Check for common non-nitrile functional groups containing C#N
    non_nitrile_groups = ["[C;$(C=O)][C;$(C#N)]", # amides
                          "[N;$(N=C)][C;$(C#N)]", # carboxamides
                          "[n;$(n:n)]=[C;$(C#N)]"] # heterocycles with nitrile substituents
    for group in non_nitrile_groups:
        pattern = Chem.MolFromSmarts(group)
        if mol.HasSubstructMatch(pattern):
            return False, f"Contains non-nitrile functional group: {group}"

    # Check molecular weight - nitriles typically < 300 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 300:
        return False, "Molecular weight too high for a simple nitrile"

    # Check atom counts - nitriles typically have < 20 non-hydrogen atoms
    non_h_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() != 1)
    if non_h_atoms > 20:
        return False, "Too many non-hydrogen atoms for a simple nitrile"

    # Check for aromatic rings
    if mol.GetNumAtoms() > 10 and mol.HasSubstructMatch(Chem.MolFromSmarts("a")):
        return True, "Contains nitrile group and aromatic ring(s)"

    # Check for aliphatic chains
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[C;$(C(-;C)(-;C)(-;C))]")):
        return True, "Contains nitrile group and aliphatic chain(s)"

    return True, "Contains nitrile group (C-C#N)"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:17116',
        'name': 'nitrile',
        'definition': 'A compound having the structure RC#N; thus a C-substituted derivative of hydrocyanic acid, HC#N. In systematic nomenclature, the suffix nitrile denotes the triply bound #N atom, not the carbon atom attached to it.',
        'parents': ['CHEBI:35190', 'CHEBI:33393']
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
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 512,
    'num_false_positives': 6,
    'num_true_negatives': 182369,
    'num_false_negatives': 5,
    'num_negatives': None,
    'precision': 0.9885057471264367,
    'recall': 0.9902785304339493,
    'f1': 0.9893917977682495,
    'accuracy': 0.9998304013042115
}